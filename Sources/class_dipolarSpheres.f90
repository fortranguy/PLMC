!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

use iso_fortran_env
use data_constants
use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use mod_physics
use class_neighbours
use class_mixingPotential
use class_spheres

implicit none

private

    type, extends(Spheres), public :: DipolarSpheres

        private
        
        ! Particles
        
        real(DP), dimension(:, :), allocatable, public :: M !< moments of all particles

        ! Potential
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :), allocatable :: Epot_real_tab !< tabulation : real short-range
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => DipolarSpheres_construct
        procedure :: destroy => DipolarSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: report => DipolarSpheres_report
        
        !> Take a snap shot of the configuration : orientations
        procedure :: snapShot_M => DipolarSpheres_snapShot_M
        
        !> Potential energy
        procedure :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure :: Epot_real_interpol => DipolarSpheres_Epot_real_interpol
        procedure :: Epot_real_pair => DipolarSpheres_Epot_real_pair
        procedure :: Epot_real => DipolarSpheres_Epot_real
        procedure :: Epot_self_delta => DipolarSpheres_Epot_self_delta
        procedure :: Epot_self => DipolarSpheres_Epot_self
        procedure :: Epot_neigh => DipolarSpheres_Epot_neigh
        procedure :: Epot_conf => DipolarSpheres_Epot_conf
        procedure :: consistTest => DipolarSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => DipolarSpheres_move
        procedure :: widom => DipolarSpheres_widom
        
    end type DipolarSpheres
    
contains

    subroutine DipolarSpheres_construct(this, shared_rCut)
    
        class(DipolarSpheres), intent(out) :: this
        real(DP), intent(in) :: shared_rCut
        
        this%name = "dipol"
    
        ! Particles
        this%radius = dipol_radius
        this%rMin = dipol_rMin
        this%Ncol = dipol_Ncol
        allocate(this%X(Dim, this%Ncol))
        allocate(this%M(Dim, this%Ncol))
        
        ! Monte-Carlo
        this%dx = dipol_dx
        this%dx_save = dipol_dx
        this%rejFix = dipol_rejFix
        this%Nadapt = dipol_Nadapt
        this%Nwidom = dipol_Nwidom
        
        ! Potential
        this%rCut = dipol_rCut
        this%dr = dipol_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%alpha = dipol_alpha        
        allocate(this%Epot_real_tab(this%iMin:this%iCut, 2))
        call this%Epot_real_init()
        
        ! Neighbours : same kind    
        call this%same%construct(this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%mix%construct(shared_rCut)
        call this%mix%alloc_cells()
        call this%mix%ini_cell_neighs()
    
    end subroutine DipolarSpheres_construct
    
    subroutine DipolarSpheres_destroy(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        if (allocated(this%X)) then
            deallocate(this%X)
        end if
        
        if (allocated(this%M)) then
            deallocate(this%M)
        end if
        
        if (allocated(this%Epot_real_tab)) then
            deallocate(this%Epot_real_tab)
        endif
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine DipolarSpheres_destroy
    
    !> Report
    
    subroutine DipolarSpheres_report(this, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
        write(report_unit, *) "    same_cell_coordMax(:) = ", this%same%cell_coordMax(:)
        write(report_unit, *) "    same_cell_Lsize(:) = ", this%same%cell_Lsize(:)
        write(report_unit, *) "    mix_cell_coordMax(:) = ", this%mix%cell_coordMax(:)
        write(report_unit, *) "    mix_cell_Lsize(:) = ", this%mix%cell_Lsize(:)
        
    end subroutine DipolarSpheres_report
    
    !> Configuration state : orientations
      
    subroutine DipolarSpheres_snapShot_M(this, snap_unit)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        do iCol = 1, this%Ncol
            write(snap_unit, *) this%M(:, iCol)
        end do    

    end subroutine DipolarSpheres_snapShot_M
    
    !> Potential energy
    
    subroutine DipolarSpheres_Epot_real_init(this)
    
        class(DipolarSpheres), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
        real(DP) :: alpha
        
        alpha = this%alpha
       
        ! cut
        do i = this%iMin, this%iCut
        
            r_i = real(i, DP)*this%dr
            
            this%Epot_real_tab(i, 1) = erfc(alpha*r_i)/r_i**3 + &
                                  2._DP*alpha/sqrt(PI) * exp(-alpha**2*r_i**2) / r_i**2
                                 
            this%Epot_real_tab(i, 2) = 3._DP*erfc(alpha*r_i)/r_i**5 + &
                                  2._DP*alpha/sqrt(PI) * (2_DP*alpha**2+3._DP/r_i**2) * &
                                                         exp(-alpha**2*r_i**2) / r_i**2
                                    
        end do
        
        ! shift        
        this%Epot_real_tab(:, 1) = this%Epot_real_tab(:, 1) - this%Epot_real_tab(this%iCut, 1)
        this%Epot_real_tab(:, 2) = this%Epot_real_tab(:, 2) - this%Epot_real_tab(this%iCut, 2)

    end subroutine DipolarSpheres_Epot_real_init
    
    !> Print the tabulated potential
    
    subroutine DipolarSpheres_Epot_real_print(this, Epot_unit)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            write(Epot_unit, *) r_i, this%Epot_real_tab(i, :)
        end do

    end subroutine DipolarSpheres_Epot_real_print

    function DipolarSpheres_Epot_real_interpol(this, r) result(Epot_real_interpol)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP), dimension(2) :: Epot_real_interpol
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            Epot_real_interpol(:) = this%Epot_real_tab(i, :) + (r-r_i)/this%dr * &
                               (this%Epot_real_tab(i+1, :) - this%Epot_real_tab(i, :))
           
        else
       
            Epot_real_interpol(:) = 0._DP
           
        end if
        
    end function DipolarSpheres_Epot_real_interpol
    
    function DipolarSpheres_Epot_real_pair(this, iCol, jCol, rVec, r) result(Epot_real_pair)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol, jCol
        real(DP), dimension(:), intent(in) :: rVec
        real(DP), intent(in) :: r
        real(DP) :: Epot_real_pair
        
        real(DP), dimension(2) :: Epot_coeff
        
        Epot_coeff(1) = dot_product(this%M(:, iCol), this%M(:, jCol))
        Epot_coeff(2) =-dot_product(this%M(:, iCol), rVec) * dot_product(this%M(:, jCol), rVec)
        
        Epot_real_pair = dot_product(Epot_coeff, this%Epot_real_interpol(r))
    
    end function DipolarSpheres_Epot_real_pair
    
    !> Total real energy
    
    function DipolarSpheres_Epot_real(this) result(Epot_real)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_real
        
        integer :: iCol, jCol
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij        
    
        Epot_real = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    rVec_ij = distVec(this%X(:, iCol), this%X(:, jCol))
                    r_ij = dot_product(rVec_ij, rVec_ij)
                    
                    Epot_real = Epot_real + this%Epot_real_pair(iCol, jCol, rVec_ij, r_ij)
                    
                end if
            end do
        end do
        
        Epot_real = 0.5_DP*Epot_real
    
    end function DipolarSpheres_Epot_real
    
    function DipolarSpheres_Epot_self_delta(this) result(Epot_self_delta)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_self_delta
    
    end function DipolarSpheres_Epot_self_delta
    
    !> Total self energy
    
    function DipolarSpheres_Epot_self(this) result(Epot_self)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_self
        
        real(DP) :: momentsSum
        integer :: iCol        
        
        momentsSum = 0._DP
        do iCol = 1, this%Ncol
            momentsSum = momentsSum + dot_product(this%M(:, iCol), this%M(:, iCol))
        end do
        
        Epot_self = 2._DP/3._DP * this%alpha**3/sqrt(PI) * momentsSum
        
    end function DipolarSpheres_Epot_self
    
    subroutine DipolarSpheres_Epot_neigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP), dimension(Dim) :: rVec
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
                
                    rVec = distVec(xCol(:), this%X(:, current%iCol))
                    r = dot_product(rVec, rVec)
                    
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    
                    energ = energ + this%Epot_real_pair(iCol, current%iCol, rVec, r)
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do            
            
        end do
    
    end subroutine DipolarSpheres_Epot_neigh
    
    !> Particle move
    
    subroutine DipolarSpheres_move(this, other, mix, same_Epot, mix_Epot, Nrej)
    
        class(DipolarSpheres), intent(inout) :: this
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix        
        real(DP), intent(inout) :: same_Epot, mix_Epot
        integer, intent(inout) :: Nrej
        
        logical :: overlap
        integer :: iOld
        real(DP) :: rand
        real(DP), dimension(Dim) :: xRand, xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: dEpot
        real(DP) :: same_dEpot, mix_dEpot
        real(DP) :: same_eNew, same_eOld
        real(DP) :: mix_eNew, mix_eOld
        
        call random_number(rand)
        iOld = int(rand*this%Ncol) + 1
        
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + (xRand(:)-0.5_DP)*this%dx(:)
        xNew(:) = modulo(xNew(:), Lsize(:))
        same_iCellNew = this%same%position_to_cell(xNew)
        call this%Epot_neigh(iOld, xNew, same_iCellNew, overlap, same_eNew)
        
        if (.not. overlap) then
        
            mix_iCellNew = this%mix%position_to_cell(xNew)
            call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%X, overlap, mix_eNew)
                        
            if (.not. overlap) then
    
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%Epot_neigh(iOld, this%X(:, iOld), same_iCellOld, overlap, same_eOld)                    
                same_dEpot = same_eNew - same_eOld
                    
                mix_iCellOld = this%mix%position_to_cell(this%X(:, iOld))
                call mix%Epot_neigh(this%X(:, iOld), mix_iCellOld, this%mix, other%X, overlap, mix_eOld)
                mix_dEpot = mix_eNew - mix_eOld
                
                dEpot = same_dEpot + mix_dEpot
                
                call random_number(rand)
                if (rand < exp(-dEpot/Tstar)) then
                
                    this%X(:, iOld) = xNew(:)
                    same_Epot = same_Epot + same_dEpot
                    mix_Epot = mix_Epot + mix_dEpot
                    
                    if (same_iCellOld /= same_iCellNew) then
                        call this%same%remove_cell_col(iOld, same_iCellOld)
                        call this%same%add_cell_col(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then
                        call other%mix%remove_cell_col(iOld, mix_iCellOld)
                        call other%mix%add_cell_col(iOld, mix_iCellNew)
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
    
    end subroutine DipolarSpheres_move
    
    !> Widom's method

    subroutine DipolarSpheres_widom(this, activ)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(inOut) :: activ 
        
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
            call this%Epot_neigh(0, xTest, iCellTest, overlap, enTest) 
            
            if (.not. overlap) then
                widTestSum = widTestSum + exp(-enTest/Tstar)
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine DipolarSpheres_widom

    !> Total potential energy
    
    function DipolarSpheres_Epot_conf(this) result(Epot_conf)
    
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
        
        Epot_conf = this%Epot_real() - this%Epot_self()
    
    end function DipolarSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine DipolarSpheres_consistTest(this, Epot, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
        real(DP) :: difference
        
        Epot_conf = this%Epot_conf()
        difference = abs((Epot_conf-Epot)/Epot_conf)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", difference
        
        if (difference > consist_tiny) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine DipolarSpheres_consistTest

end module class_dipolarSpheres
