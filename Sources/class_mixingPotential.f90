!> \brief Description of the  Mixing Potential class

module class_mixingPotential

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_box, only: Ndim
use data_particles, only: mix_delta
use data_potential, only: mix_rMin_factor, mix_rCut, mix_dr, mix_epsilon, mix_alpha
use data_neighbourCells, only: NnearCell
use module_types, only: Node
use module_physics_micro, only: set_discrete_length, dist_PBC, Epot_yukawa
use class_neighbourCells
use class_hardSpheres

implicit none

private

    type, public :: MixingPotential

        private
        
        character(len=5) :: name
        
        real(DP) :: delta !< demixing length
        real(DP) :: sigma
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: rCut !< short-range cut
        real(DP) :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation: minimum distance
        integer :: iCut !< maximum index of tabulation: until potential cut
        real(DP) :: epsilon !< factor in Yukawa
        real(DP) :: alpha !< coefficient in Yukawa
        real(DP), dimension(:), allocatable :: Epot_tab !< tabulation
        
        real(DP), dimension(Ndim) :: cell_size

    contains
    
        procedure :: construct => MixingPotential_construct
        procedure :: destroy => MixingPotential_destroy
        procedure :: print_report => MixingPotential_print_report

        procedure :: get_sigma => MixingPotential_get_sigma
        procedure :: get_rCut => MixingPotential_get_rCut
        procedure :: set_cell_size => MixingPotential_set_cell_size
        procedure :: get_cell_size => MixingPotential_get_cell_size
        
        procedure :: test_overlap => MixingPotential_test_overlap
        procedure, private :: Epot_set_tab => MixingPotential_Epot_set_tab
        procedure :: Epot_set => MixingPotential_Epot_set
        procedure :: Epot_print => MixingPotential_Epot_print
        procedure, private :: Epot_pair => MixingPotential_Epot_pair
        procedure :: Epot_neighCells => MixingPotential_Epot_neighCells
        procedure :: Epot_conf => MixingPotential_Epot_conf

    end type

contains

    subroutine MixingPotential_construct(this, type1_sigma, type2_sigma)

        class(MixingPotential), intent(out) :: this
        real(DP), intent(in) :: type1_sigma, type2_sigma
        
        this%name = "[mix]"
        write(output_unit, *) this%name, " class construction"
        
        ! Particles
        this%delta = mix_delta
        this%sigma = (type1_sigma + type2_sigma)/2._DP + this%delta

    end subroutine MixingPotential_construct

    subroutine MixingPotential_destroy(this)
    
        class(MixingPotential), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%Epot_tab)) deallocate(this%Epot_tab)
    
    end subroutine MixingPotential_destroy
    
    !> Report
    
    subroutine MixingPotential_print_report(this, report_unit)
    
        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit, *) "    delta = ", this%delta
        
        write(report_unit, *) "    epsilon = ", this%epsilon
        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
    end subroutine MixingPotential_print_report

    !> Accessors & Mutators

    pure function MixingPotential_get_sigma(this) result(get_sigma)
        class(MixingPotential), intent(in) :: this
        real(DP) :: get_sigma
        get_sigma = this%sigma
    end function MixingPotential_get_sigma
    
    pure function MixingPotential_get_rCut(this) result(get_rCut)
        class(MixingPotential), intent(in) :: this
        real(DP) :: get_rCut
        get_rCut = this%rCut
    end function MixingPotential_get_rCut
    
    pure subroutine MixingPotential_set_cell_size(this)
        class(MixingPotential), intent(inout) :: this
        this%cell_size(:) = this%rCut
    end subroutine MixingPotential_set_cell_size
    
    pure function MixingPotential_get_cell_size(this) result(get_cell_size)
        class(MixingPotential), intent(in) :: this
        real(DP), dimension(Ndim) :: get_cell_size
        get_cell_size(:) = this%cell_size(:)
   end function MixingPotential_get_cell_size
    
    !> Overlapt test
    
    subroutine MixingPotential_test_overlap(this, type1, type2)
    
        class(MixingPotential), intent(in) :: this
        class(HardSpheres), intent(in) :: type1, type2
        
        integer :: type1_iCol, type2_iCol
        real(DP) :: r_mix
        
        do type1_iCol = 1, type1%get_Ncol()
            do type2_iCol = 1, type2%get_Ncol()
                    
                r_mix = dist_PBC(type1%positions(:, type1_iCol), type2%positions(:, type2_iCol))
                if (r_mix < this%rMin) then
                    write(error_unit, *) this%name, ":    Overlap !", type1_iCol, type2_iCol
                    write(error_unit, *) "    r_mix = ", r_mix
                    error stop
                end if

            end do
        end do

        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine MixingPotential_test_overlap

    !> MixingPotential energy
    !> Tabulation of Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_{min})}}{r} \f]
    
    pure subroutine MixingPotential_Epot_set_tab(this)
    
        class(MixingPotential), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
       
        ! cut
        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            this%Epot_tab(i) = Epot_yukawa(this%epsilon, this%alpha, this%rMin, r_i)
        end do
        
        ! shift
        this%Epot_tab(:) = this%Epot_tab(:) - this%Epot_tab(this%iCut)

    end subroutine MixingPotential_Epot_set_tab
    
    subroutine MixingPotential_Epot_set(this)
    
        class(MixingPotential), intent(inout) :: this

        this%rMin = mix_rMin_factor * this%sigma
        this%rCut = mix_rCut
        
        if (this%rCut < this%rMin) then
            this%rCut = this%rMin
        end if
        
        this%dr = mix_dr
        call set_discrete_length(this%rMin, this%dr)
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr) + 1
        this%epsilon = mix_epsilon
        this%alpha = mix_alpha
        
        if (allocated(this%Epot_tab)) deallocate(this%Epot_tab)
        allocate(this%Epot_tab(this%iMin:this%iCut))
        call this%Epot_set_tab()
        
    end subroutine MixingPotential_Epot_set
    
    !> Print the tabulated potential
    
    subroutine MixingPotential_Epot_print(this, Epot_unit)

        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            write(Epot_unit, *) r_i, this%Epot_tab(i)
        end do

    end subroutine MixingPotential_Epot_print

    pure function MixingPotential_Epot_pair(this, r) result(Epot_pair)
        
        class(MixingPotential), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP) :: Epot_pair
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%rCut) then
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            Epot_pair = this%Epot_tab(i) + (r-r_i)/this%dr * (this%Epot_tab(i+1)-this%Epot_tab(i))
        else
            Epot_pair = 0._DP
        end if
        
    end function MixingPotential_Epot_pair
    
    subroutine MixingPotential_Epot_neighCells(this, iCol, xCol, iTotalCell, neighCells, &
                                               other_positions, overlap, energ)
        
        class(MixingPotential), intent(in) :: this
        integer, intent(in) :: iCol
        real(DP), dimension(:), intent(in) :: xCol !< type A
        integer, intent(in) :: iTotalCell !< type A in mix grid
        type(NeighbourCells), intent(in) :: neighCells
        real(DP), dimension(:, :), contiguous, intent(in) :: other_positions
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNearCell = 1, NnearCell
        
            nearCell_index = neighCells%near_among_total(iNearCell, iTotalCell)
            current => neighCells%beginCells(nearCell_index)%particle%next
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next

                if (current%iCol /= iCol) then
                    r = dist_PBC(xCol(:), other_positions(:, current%iCol))
                    if (r < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    energ = energ + this%Epot_pair(r)
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine MixingPotential_Epot_neighCells
    
    !> Total potential energy
    
    pure function MixingPotential_Epot_conf(this, type1, type2) result(Epot_conf)
    
        class(MixingPotential), intent(in) :: this
        class(HardSpheres), intent(in) :: type1, type2
        real(DP) :: Epot_conf

        integer :: type1_iCol, type2_iCol
        real(DP) :: r_mix

        Epot_conf = 0._DP
        
        if (this%epsilon < real_zero) then
            return
        end if
        
        do type1_iCol = 1, type1%get_Ncol()
            do type2_iCol = 1, type2%get_Ncol()
                
                r_mix = dist_PBC(type1%positions(:, type1_iCol), type2%positions(:, type2_iCol))
                Epot_conf = Epot_conf + this%Epot_pair(r_mix)

            end do
        end do
    
    end function MixingPotential_Epot_conf

end module class_mixingPotential
