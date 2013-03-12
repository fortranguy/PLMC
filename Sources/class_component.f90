module mod_pbc

use data_constants
use data_cell

implicit none

    contains

! Distance between 2 positions : PBC ------------------------------------------
    
    function dist(X1, X2)
    
        real(DP), dimension(Dim), intent(in) :: X1, X2
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: dist
        
        DeltaX(:) = X2(:) - X1(:)
        DeltaX(:) = modulo(DeltaX(:), Lsize(:))
        
        where( DeltaX(:) > LsizeMi(:) )
            DeltaX(:) = DeltaX(:) - Lsize(:)
        end where
        
        dist = sqrt(dot_product(DeltaX, DeltaX))
    
    end function dist

end module mod_pbc

!***********************************************************************
!* MODULE: Component class                                             *
!***********************************************************************

module class_component

use data_cell
use data_particles
use data_potentiel
use data_neighbours
use data_mc
use mod_pbc
use class_neighbours

implicit none

private
public :: sph_constructor

    type, public :: Component

        ! Particles

        real(DP), private :: radius
        real(DP), private :: rmin
        integer, private ::  Ncol !!!
        real(DP), dimension(:, :), allocatable :: X !!!

        ! Monte-Carlo
        
        real(DP), dimension(Dim), private :: dx

        ! Potential

        real(DP), private :: rcut !!!
        real(DP), private :: pas
        integer, private :: iMin
        integer, private :: Ntab
        real(DP), private :: epsilon
        real(DP), private :: alpha
        real(DP), dimension(:), allocatable, private :: Vtab
        
    contains
    
        procedure :: rapport => component_rapport
        procedure :: destructor => component_destructor
        procedure :: snapShot => component_snapShot
        procedure :: overlapTest => component_overlapTest
        
        procedure :: adapt_dx => component_adapt_dx
        procedure :: getDx => component_getDx
        
        procedure :: ePotIni => component_ePotIni
        procedure :: ePot => component_ePot
        procedure :: ePotNeigh => component_ePotNeigh
        procedure :: enTotCalc => component_enTotCalc
        
        procedure :: mcMove => component_mcMove
        procedure :: widom => component_widom
        
    end type Component
    
contains

    function sph_constructor()
    
        type(Component) :: sph_constructor
    
        ! Construction                

        sph_constructor%radius = sph_radius
        sph_constructor%rmin = sph_rmin
        sph_constructor%Ncol = sph_Ncol
        allocate(sph_constructor%X(Dim, sph_Ncol))
        
        sph_constructor%dx = sph_dx
        
        sph_constructor%rcut = sph_rcut
        sph_constructor%pas = sph_pas
        sph_constructor%iMin = sph_iMin
        sph_constructor%Ntab = sph_Ntab        
        sph_constructor%epsilon = sph_epsilon
        sph_constructor%alpha = sph_alpha        
        allocate(sph_constructor%Vtab(sph_iMin:sph_Ntab))
        call sph_constructor%ePotIni()
        
        sph_constructor%cell_Lsize(:) = [sph_rcut, sph_rcut, sph_rcut]
        sph_constructor%cell_coordMax(:) = int(Lsize(:)/sph_rcut)
        allocate(sph_constructor%cell_neighs(cell_neighs_nb, &
            product( int(Lsize(:)/sph_rcut) )))
    
    end function sph_constructor
    
    ! Report ------------------------------------------------------------------
    
    subroutine component_rapport(this, nWidom, unitRapport)
    
        class(Component), intent(in) :: this
        integer, intent(in) :: nWidom
        integer, intent(in) :: unitRapport    
        
        write(unitRapport, *) "Simulation MC_C :"
        write(unitRapport ,*) "    Lsize(:) = ", Lsize(:)
        write(unitRapport ,*) "    Vol = ", product(Lsize)
        write(unitRapport ,*) "    Ncol = ", this%Ncol
        write(unitRapport ,*) "    nWidom = ", nWidom
        write(unitRapport, *) "    Nstep = ", Nstep
        write(unitRapport, *) "    Ntherm = ", Ntherm
        write(unitRapport, *) "    Nmove = ", Nmove
        write(unitRapport, *) "    epsilon = ", this%epsilon
        write(unitRapport, *) "    alpha = ", this%alpha
        write(unitRapport, *) "    rcut = ", this%rcut
        write(unitRapport, *) "    pas = ", this%pas
        write(unitRapport, *) "    cell_coordMax(:) = ", this%cell_coordMax(:)
        write(unitRapport, *) "    cell_Lsize(:) = ", this%cell_Lsize(:)
        
    end subroutine component_rapport
    
    subroutine component_destructor(this)
    
        class(Component), intent(inout) :: this
        
        deallocate(this%X)
        deallocate(this%Vtab)
    
    end subroutine component_destructor
    
    ! Configuration state -----------------------------------------------------
      
    subroutine component_snapShot(this, unitSnap)
        
        class(Component), intent(in) :: this
        integer, intent(in) :: unitSnap
    
        integer :: iCol
        
        do iCol = 1, this%Ncol
            write(unitSnap, *) this%X(:, iCol)
        end do    

    end subroutine component_snapShot
    
    ! Overlapt test -----------------------------------------------------------
    
    subroutine component_overlapTest(this)
    
        class(Component), intent(in) :: this
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist(this%X(:, iCol), this%X(:, jCol))
                    if (r_ij < this%rmin) then
                        write(*, *) "    Overlap !", iCol, jCol
                        write(*, * ) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(*, *) "    Overlap test : OK !"
    
    end subroutine component_overlapTest
    
    ! Adaptation of dx during the thermalisation ------------------------------
    
    subroutine component_adapt_dx(this, iStep, tauxRejectsSum, unitRapport)
    
         class(Component), intent(inout) :: this 
        integer, intent(in) :: iStep, unitRapport
        real(DP), intent(in) :: tauxRejectsSum    
        
        integer, parameter :: multiple = 2**2
        real(DP) :: tauxRejects
        real(DP), parameter :: tauxRejectsFix = 0.5_DP
        real(DP), parameter :: dx_eps = 0.05_DP, taux_eps = 0.05_DP
        real(DP), parameter :: more = 1._DP+dx_eps, less = 1._DP-dx_eps
        
        tauxRejects = 0._DP
        
        if (mod(iStep, multiple) == 0 .and. iStep>2) then
        
            tauxRejects = tauxRejectsSum/real(iStep-1, DP)
        
            if (tauxRejects < tauxRejectsFix - taux_eps) then            
                this%dx(:) = this%dx(:) * more
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            else if (tauxRejects > tauxRejectsFix + taux_eps) then
                this%dx(:) = this%dx(:) * less
                this%dx(:) = modulo(this%dx(:), Lsize(:))
            end if

        end if
        
        if (iStep == Ntherm) then
        
            if (tauxRejects == 0._DP) then
                write(*, *) "Problème adaptation dx."
                stop
            end if
            
            write(unitRapport, *) "Déplacement :"
            write(unitRapport, *) "    dx(:) = ", this%dx(:)
            write(unitRapport, *) "    écart relatif rejet = ", &
                abs(tauxRejects - tauxRejectsFix)/tauxRejectsFix
            
        end if
    
    end subroutine component_adapt_dx
    
    ! -----------------------
    
    function component_getDx(this)
        
        class(Component), intent(in) :: this
        
        real(DP) :: component_getDx
        
        component_getDx = sqrt(dot_product(this%dx, this%dx))
        
    end function component_getDx
    
    ! Potential energy --------------------------------------------------------
    
    subroutine component_ePotIni(this)
    
        class(Component), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%Ntab       
            r_i = real(i, DP)*this%pas
            this%Vtab(i) = this%epsilon * exp(-this%alpha*(r_i-this%rmin))/r_i           
        end do
        
        ! shift        
        this%Vtab(:) = this%Vtab(:) - this%epsilon * &
            exp(-this%alpha*(this%rcut-this%rmin)) / this%rcut

    end subroutine component_ePotIni

    function component_ePot(this, r) result(ePot)
        
        class(Component), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i, ePot
       
        if (r < this%rcut) then
       
            i = int(r/this%pas)
            r_i = real(i, DP)*this%pas
            ePot = this%Vtab(i) + (r-r_i)/this%pas * &
                (this%Vtab(i+1)-this%Vtab(i))
           
        else
       
            ePot = 0.
           
        end if
        
    end function component_ePot
    
    ! -----------------------
    
    subroutine component_ePotNeigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(Component), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP) :: r
    
        type(Link), pointer :: courant => null(), suivant => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = this%cell_neighs(iNeigh, iCell)
            courant => this%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(courant%next)) cycle
            
            do
            
                suivant => courant%next
            
                if (courant%iCol /= iCol) then
                
                    r = dist(xCol(:), this%X(:, courant%iCol))
                    if (r < this%rmin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%ePot(r)
       
                end if
                
                if (.not. associated(suivant%next)) exit
                
                courant => suivant
            
            end do            
            
        end do
    
    end subroutine component_ePotNeigh
    
    ! Particle move -----------------------------------------------------------
    
    subroutine component_mcMove(this, enTot, Nrejects)
    
        class(Component), intent(inout) :: this
        real(DP), intent(inout) :: enTot
        integer, intent(inout) :: Nrejects
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xNew
        integer :: iCellBefore, iCellAfter
        real(DP) :: eNew, eOld, dEn
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xNew)
        xNew(:) = this%X(:, iOld) + (xNew(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        iCellAfter = this%position_to_cell(xNew)
        call this%ePotNeigh(iOld, xNew, iCellAfter, overlap, eNew)
        
        if (.not. overlap) then
        
            iCellBefore = this%position_to_cell(this%X(:, iOld))
            call this%ePotNeigh(iOld, this%X(:, iOld), iCellBefore, overlap, &
                eOld)
            
            dEn = eNew - eOld
        
            call random_number(rand)
            if ( rand < exp(-dEn/Tstar) ) then
                this%X(:, iOld) = xNew(:)
                enTot = enTot + dEn
                
                if ( iCellBefore /= iCellAfter ) then                
                    call this%remove_cell_col(iOld, iCellBefore)
                    call this%add_cell_col(iOld, iCellAfter)
                end if
                
            else
                Nrejects = Nrejects + 1
            end if
            
        else
        
            Nrejects = Nrejects + 1
            
        end if
    
    end subroutine component_mcMove
    
    ! Widom's method -----------------------------------------------------------

    subroutine component_widom(this, nWidom, activExInv)
        
        class(Component), intent(in) :: this
        integer, intent(in) :: nWidom
        real(DP), intent(inOut) :: activExInv 
        
        integer :: iWid
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest
        integer :: iCellTest
        logical :: overlap        
        real(DP) :: enTest
        
        widTestSum = 0._DP
        
        do iWid = 1, nWidom           
            
            call random_number(xTest)
            xTest(:) = Lsize(:) * xTest(:)    
            iCellTest = this%position_to_cell(xTest)
            call this%ePotNeigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activExInv = widTestSum/real(nWidom, DP)
        
    end subroutine component_widom

    ! Total potential energy
    
    function component_enTotCalc(this) result(enTotCalc)
    
        class(Component), intent(in) :: this
        
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
    
    end function component_enTotCalc

end module class_component
