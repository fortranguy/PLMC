!> \brief Description of the Interacting Sphere class

module class_interacting

use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use class_neighbours
use mod_physics
use class_spheres

implicit none

private
public :: inter_constructor

    type, extends(Spheres), public :: Interacting

        private

        ! Potential :
        real(DP)  :: pas !< discretisation step
        integer :: iMin !< minimum index of tabulation
        integer :: Ntab !< maximum index of tabulation
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: Vtab !< tabulation
        
    contains

        !> Destructor of the class
        procedure :: destructor => Interacting_destructor
        
        !> Print a report of the component in a file
        procedure :: report => Interacting_report
              
        !> Potential energy
        procedure :: ePotIni => Interacting_ePotIni
        procedure :: ePot => Interacting_ePot
        procedure :: ePotNeigh => Interacting_ePotNeigh
        procedure :: enTotCalc => Interacting_enTotCalc
        
        !> Monte-Carlo
        procedure :: mcMove => Interacting_mcMove
        procedure :: widom => Interacting_widom
        
    end type Interacting
    
contains

    function inter_constructor()
    
        type(Interacting) :: inter_constructor
    
        ! Construction                

        inter_constructor%radius = inter_radius
        inter_constructor%rMin = inter_rMin
        inter_constructor%Ncol = inter_Ncol
        allocate(inter_constructor%X(Dim, inter_Ncol))
        
        inter_constructor%dx = inter_dx
        
        inter_constructor%rCut = inter_rCut
        inter_constructor%pas = inter_pas
        inter_constructor%iMin = inter_iMin
        inter_constructor%Ntab = inter_Ntab        
        inter_constructor%epsilon = inter_epsilon
        inter_constructor%alpha = inter_alpha        
        allocate(inter_constructor%Vtab(inter_iMin:inter_Ntab))
        call inter_constructor%ePotIni()
        
        !	Neighbours        
        inter_constructor%same = neigh_constructor(inter_rCut)
        call inter_constructor%same%alloc_cells()
        call inter_constructor%same%ini_cell_neighs()
    
    end function inter_constructor
    
    subroutine Interacting_destructor(this)
    
        class(Interacting), intent(inout) :: this
        
        deallocate(this%X)
        deallocate(this%Vtab)
        call this%same%destructor()
    
    end subroutine Interacting_destructor
    
    !> Report
    
    subroutine Interacting_report(this, nWidom, unitReport)
    
        class(Interacting), intent(in) :: this
        integer, intent(in) :: nWidom
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitReport ,*) "    Vol = ", product(Lsize)
        write(unitReport ,*) "    Ncol = ", this%Ncol
        write(unitReport ,*) "    nWidom = ", nWidom
        write(unitReport, *) "    Nstep = ", Nstep
        write(unitReport, *) "    Ntherm = ", Ntherm
        write(unitReport, *) "    Nmove = ", Nmove
        write(unitReport, *) "    epsilon = ", this%epsilon
        write(unitReport, *) "    alpha = ", this%alpha
        write(unitReport, *) "    rCut = ", this%rCut
        write(unitReport, *) "    pas = ", this%pas
        write(unitReport, *) "    cell_coordMax(:) = ", &
        	this%same%cell_coordMax(:)
        write(unitReport, *) "    cell_Lsize(:) = ", this%same%cell_Lsize(:)
        
    end subroutine Interacting_report
    
    !> Potential energy
    !> Tabulation of Yukawa potential    
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine Interacting_ePotIni(this)
    
        class(Interacting), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%Ntab       
            r_i = real(i, DP)*this%pas
            this%Vtab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rMin))/r_i           
        end do
        
        ! shift        
        this%Vtab(:) = this%Vtab(:) - this%epsilon * &
            exp(-this%alpha*(this%rCut-this%rMin)) / this%rCut

    end subroutine Interacting_ePotIni

    function Interacting_ePot(this, r) result(ePot)
        
        class(Interacting), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < this%rCut) then
       
            i = int(r/this%pas)
            r_i = real(i, DP)*this%pas
            ePot = this%Vtab(i) + (r-r_i)/this%pas * &
                (this%Vtab(i+1)-this%Vtab(i))
           
        else
       
            ePot = 0._DP
           
        end if
        
    end function Interacting_ePot
    
    subroutine Interacting_ePotNeigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(Interacting), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = this%same%cell_neighs(iNeigh, iCell)
            current => this%same%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= iCol) then
                
                    r = dist(xCol(:), this%X(:, current%iCol))
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%ePot(r)
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine Interacting_ePotNeigh
    
    !> Particle move
    
    subroutine Interacting_mcMove(this, enTot, Nrejects)
    
        class(Interacting), intent(inout) :: this
        real(DP), intent(inout) :: enTot
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: iCellBefore, iCellAfter
        real(DP) :: eNew, eOld, dEn
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        iCellAfter = this%same%position_to_cell(xNew)
        call this%ePotNeigh(iOld, xNew, iCellAfter, overlap, eNew)
        
        if (.not. overlap) then
        
            iCellBefore = this%same%position_to_cell(this%X(:, iOld))
            call this%ePotNeigh(iOld, this%X(:, iOld), iCellBefore, overlap, &
                eOld)
            
            dEn = eNew - eOld
        
            call random_number(rand)
            if ( rand < exp(-dEn/Tstar) ) then
                this%X(:, iOld) = xNew(:)
                enTot = enTot + dEn
                
                if ( iCellBefore /= iCellAfter ) then                
                    call this%same%remove_cell_col(iOld, iCellBefore)
                    call this%same%add_cell_col(iOld, iCellAfter)
                end if
                
            else
                Nrejects = Nrejects + 1
            end if
            
        else
        
            Nrejects = Nrejects + 1
            
        end if
    
    end subroutine Interacting_mcMove
    
    !> Widom's method

    subroutine Interacting_widom(this, nWidom, activExInv)
        
        class(Interacting), intent(in) :: this
        integer, intent(in) :: nWidom
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            iCellTest = this%same%position_to_cell(xTest)
            call this%ePotNeigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine Interacting_widom

    !> Total potential energy
    
    function Interacting_enTotCalc(this) result(enTotCalc)
    
        class(Interacting), intent(in) :: this
        
        integer :: iCol, jCol
        real(DP) :: r_ij
        real(DP) :: enTotCalc
    
        enTotCalc = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                
                    r_ij = dist(this%X(:, iCol), this%X(:, jCol))
                    enTotCalc = enTotCalc + this%ePot(r_ij)
                    
                end if
            end do
        end do
        
        enTotCalc = 0.5_DP*enTotCalc
    
    end function Interacting_enTotCalc

end module class_interacting
