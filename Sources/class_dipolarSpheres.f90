!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

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
        
        real(DP), dimension(:, :), allocatable :: M !< moments of all particles

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
        
        !> Potential energy
        procedure :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure :: Epot_real_pair => DipolarSpheres_Epot_real_pair
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

    function DipolarSpheres_Epot_real_pair(this, r) result(Epot_real_pair)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        
        integer :: i
        real(DP) :: r_i
        real(DP), dimension(2) :: Epot_real_pair
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            Epot_real_pair(:) = this%Epot_real_tab(i, :) + &
                           (r-r_i)/this%dr * (this%Epot_real_tab(i+1, :) - this%Epot_real_tab(i, :))
           
        else
       
            Epot_real_pair(:) = 0._DP
           
        end if
        
    end function DipolarSpheres_Epot_real_pair
    
    subroutine DipolarSpheres_Epot_neigh(this, iCol, xCol, iCell, overlap, energ)
        
        class(DipolarSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(Dim), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP), dimension(Dim) :: r_vec
        real(DP) :: r
        real(DP) :: Epot_real
        real(DP), dimension(2) :: Epot_real_coeff
        
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
                
                    r_vec = dist_vec(xCol(:), this%X(:, current%iCol))
                    r = dot_product(r_vec, r_vec)
                    
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    
                    Epot_real_coeff(1) = dot_product(this%M(:, iCol), this%M(:, current%iCol))
                    Epot_real_coeff(2) =-dot_product(this%M(:, iCol), r_vec) * &
                                         dot_product(this%M(:, current%iCol), r_vec)
                    
                    Epot_real = dot_product(Epot_real_coeff, this%Epot_real_pair(r))
                     
                    energ = energ + Epot_real
       
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
        
        integer :: iCol, jCol
        real(DP), dimension(Dim) :: r_vec
        real(DP) :: r_ij
        real(DP) :: Epot_real
        real(DP), dimension(2) :: Epot_real_coeff
    
        Epot_conf = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    r_vec = dist_vec(this%X(:, iCol), this%X(:, jCol))
                    r_ij = dot_product(r_vec, r_vec)
                    
                    Epot_real_coeff(1) = dot_product(this%M(:, iCol), this%M(:, jCol))
                    Epot_real_coeff(2) =-dot_product(this%M(:, iCol), r_vec) * &
                                         dot_product(this%M(:, jCol), r_vec)
                    
                    Epot_real = dot_product(Epot_real_coeff, this%Epot_real_pair(r_ij))
                    
                    Epot_conf = Epot_conf + Epot_real
                    
                end if
            end do
        end do
        
        Epot_conf = 0.5_DP*Epot_conf
    
    end function DipolarSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine DipolarSpheres_consistTest(this, Epot, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
    
        Epot_conf = this%Epot_conf()
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", abs((Epot_conf-Epot)/Epot_conf)
    
    end subroutine DipolarSpheres_consistTest

end module class_dipolarSpheres
