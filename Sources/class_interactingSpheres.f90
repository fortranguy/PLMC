!> \brief Description of the InteractingSpheres class

module class_interactingSpheres

use data_constants
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

    type, extends(Spheres), public :: InteractingSpheres

        private 

        ! Potential :
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: ePot_tab !< tabulation
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => InteractingSpheres_construct
        procedure :: destroy => InteractingSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: report => InteractingSpheres_report
              
        !> Potential energy
        procedure :: ePot_init => InteractingSpheres_ePot_init
        procedure :: ePot => InteractingSpheres_ePot
        procedure :: ePot_neigh => InteractingSpheres_ePot_neigh
        procedure :: ePot_total => InteractingSpheres_ePot_total
        procedure :: consistTest => InteractingSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => InteractingSpheres_move
        procedure :: widom => InteractingSpheres_widom
        
    end type InteractingSpheres
    
contains

    subroutine InteractingSpheres_construct(this)
    
        class(InteractingSpheres), intent(out) :: this
        
        this%name = "inter"
    
        ! Particles
        this%radius = inter_radius
        this%rMin = inter_rMin
        this%Ncol = inter_Ncol
        allocate(this%X(Dim, this%Ncol))
        
        ! Monte-Carlo
        this%dx = inter_dx
        this%dx_save = inter_dx
        this%Nadapt = inter_Nadapt
        this%Nwidom = inter_Nwidom
        
        ! Potential
        this%rCut = inter_rCut
        this%dr = inter_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%epsilon = inter_epsilon
        this%alpha = inter_alpha        
        allocate(this%ePot_tab(this%iMin:this%iCut))
        call this%ePot_init()
        
        ! Neighbours        
        call this%same%construct(inter_rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
    
    end subroutine InteractingSpheres_construct
    
    subroutine InteractingSpheres_destroy(this)
    
        class(InteractingSpheres), intent(inout) :: this
        
        deallocate(this%X)
        deallocate(this%ePot_tab)
        call this%same%destroy()
    
    end subroutine InteractingSpheres_destroy
    
    !> Report
    
    subroutine InteractingSpheres_report(this, unitReport)
    
        class(InteractingSpheres), intent(in) :: this
        integer, intent(in) :: unitReport    
        
        write(unitReport, *) "Simulation MC_C :"
        write(unitReport ,*) "    Ncol = ", this%Ncol
        write(unitReport ,*) "    Nwidom = ", this%Nwidom
        write(unitReport, *) "    epsilon = ", this%epsilon
        write(unitReport, *) "    alpha = ", this%alpha
        write(unitReport, *) "    rCut = ", this%rCut
        write(unitReport, *) "    dr = ", this%dr
        write(unitReport, *) "    cell_coordMax(:) = ", &
        	this%same%cell_coordMax(:)
        write(unitReport, *) "    cell_Lsize(:) = ", this%same%cell_Lsize(:)
        
    end subroutine InteractingSpheres_report
    
    !> Potential energy
    !> Tabulation of Yukawa potential    
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    subroutine InteractingSpheres_ePot_init(this)
    
        class(InteractingSpheres), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut       
            r_i = real(i, DP)*this%dr
            this%ePot_tab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rMin))&
            /r_i
        end do
        
        ! shift        
        this%ePot_tab(:) = this%ePot_tab(:) - this%epsilon * &
            exp(-this%alpha*(this%rCut-this%rMin)) / this%rCut

    end subroutine InteractingSpheres_ePot_init

    function InteractingSpheres_ePot(this, r) result(ePot)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            ePot = this%ePot_tab(i) + (r-r_i)/this%dr * &
                (this%ePot_tab(i+1)-this%ePot_tab(i))
           
        else
       
            ePot = 0._DP
           
        end if
        
    end function InteractingSpheres_ePot
    
    subroutine InteractingSpheres_ePot_neigh(this, iCol, xCol, iCell, overlap, &
    	energ)
        
        class(InteractingSpheres), intent(in) :: this        
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
    
    end subroutine InteractingSpheres_ePot_neigh
    
    !> Particle move
    
    subroutine InteractingSpheres_move(this, ePot_total, Nrej)
    
        class(InteractingSpheres), intent(inout) :: this
        real(DP), intent(inout) :: ePot_total
        integer, intent(inout) :: Nrej
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: other_iCellNew
        real(DP) :: eNew, eOld, dEn
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%ePot_neigh(iOld, xNew, same_iCellNew, overlap, eNew)
        
        if (.not. overlap) then
        
            other_iCellNew = this%other%position_to_cell(xNew)
            call mix%ePot_neigh(xNew, other_iCellNew, this%other, other_X, &
                overlap, other_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%ePot_neigh(iOld, this%X(:, iOld), same_iCellOld, &
                    overlap, eOld)
                    
                other_iCellOld = this%other%position_to_cell(this%X(:, iOld))
                call mix%ePot_neigh(this%X(:, iOld), other_iCellNew, &
                    overlap, other_eOld)
                
                dEn = eNew - eOld
            
                call random_number(rand)
                if ( rand < exp(-dEn/Tstar) ) then
                    this%X(:, iOld) = xNew(:)
                    ePot_total = ePot_total + dEn
                    
                    if ( same_iCellOld /= same_iCellNew ) then                
                        call this%same%remove_cell_col(iOld, same_iCellOld)
                        call this%same%add_cell_col(iOld, same_iCellNew)
                    end if
                    
                else
                    Nrej = Nrej + 1
                end if
         
            else
                Nrej = Nrej + 1                
            end if
            
            
        else
        
            Nrej = Nrej + 1
            
        end if
    
    end subroutine InteractingSpheres_move
    
    !> Widom's method

    subroutine InteractingSpheres_widom(this, activExInv)
        
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xRand, xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)    
            iCellTest = this%same%position_to_cell(xTest)
            call this%ePot_neigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(this%Nwidom, DP)
        
    end subroutine InteractingSpheres_widom

    !> Total potential energy
    
    function InteractingSpheres_ePot_total(this) result(ePot_total)
    
        class(InteractingSpheres), intent(in) :: this
        
        integer :: iCol, jCol
        real(DP) :: r_ij
        real(DP) :: ePot_total
    
        ePot_total = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                
                    r_ij = dist(this%X(:, iCol), this%X(:, jCol))
                    ePot_total = ePot_total + this%ePot(r_ij)
                    
                end if
            end do
        end do
        
        ePot_total = 0.5_DP*ePot_total
    
    end function InteractingSpheres_ePot_total
    
    !> Consistency test 
    
    subroutine InteractingSpheres_consistTest(this, ePot_total_mc, unitReport)
    
        class(InteractingSpheres), intent(in) :: this
        real(DP), intent(in) :: ePot_total_mc
        integer, intent(in) :: unitReport
        
        real(DP) :: ePot_total
    
        ePot_total = this%ePot_total()
        write(unitReport, *) "Consistency test:"
        write(unitReport, *) "    ePot_total_mc = ", ePot_total_mc
        write(unitReport, *) "    ePot_total_final = ", ePot_total
        write(unitReport, *) "    relative difference = ", &
            abs(ePot_total-ePot_total_mc)/ePot_total
    
    end subroutine InteractingSpheres_consistTest

end module class_interactingSpheres
