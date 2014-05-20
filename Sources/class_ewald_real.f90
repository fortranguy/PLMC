module class_ewald_real

use data_precisions, only: DP
use data_constants, only: PI
use json_module, only: json_file
use module_physics_micro, only: set_discrete_length, ewald_real_B, ewald_real_C
use module_data, only: test_data_found

implicit none

private

    type, public :: Ewald_Real
    
        real(DP) :: min_distance
        real(DP) :: range_cut
        real(DP) :: delta
        
        integer :: i_min_distance
        integer :: i_range_cut        
        real(DP), dimension(:, :), allocatable :: tabulation
    
    contains
    
        procedure :: construct => Ewald_Real_construct
        procedure, private :: set_parameters => Ewald_Real_set_parameters
        procedure, private :: set_tabulation => Ewald_Real_set_tabulation
        procedure :: destroy => Ewald_Real_destroy
        !procedure, private :: set => Ewald_Real_set
        !procedure :: write => Ewald_Real_write
        !procedure, private :: interpolation => Ewald_Real_interpolation
        !procedure, private :: pair => Ewald_Real_pair
        !procedure :: solo => Ewald_Real_solo
        !procedure, private :: Epot_real => Ewald_Real_Epot_real
    
    end type Ewald_Real

contains

    subroutine Ewald_Real_construct(this, Box_size, alpha, min_distance, json)
        class(Ewald_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: json

        call this%set_parameters(Box_size, min_distance, json)

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        allocate(this%tabulation(this%i_min_distance:this%i_range_cut, 2))

        call this%set_tabulation(alpha)

    end subroutine Ewald_Real_construct

    !> Initialisation

    subroutine Ewald_Real_set_parameters(this, Box_size, min_distance, json)
        class(Ewald_Real), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: min_distance
        type(json_file), intent(inout) :: json

        character(len=4096) :: data_name
        logical :: found

        real(DP) :: range_cut_factor
        
        this%min_distance = min_distance

        data_name = "Potential.Dipoles.Ewald summation.real.range cut factor"
        call json%get(data_name, range_cut_factor, found)
        call test_data_found(data_name, found)
        this%range_cut = range_cut_factor * Box_size(1)

        data_name = "Potential.Dipoles.Ewald summation.real.delta"
        call json%get(data_name, this%delta, found)
        call test_data_found(data_name, found)
        call set_discrete_length(this%min_distance, this%delta)
        
        this%i_min_distance = int(this%min_distance/this%delta)
        this%i_range_cut = int(this%range_cut/this%delta) + 1

    end subroutine Ewald_Real_set_parameters

    !> Initialisation: look-up (tabulation) table

    pure subroutine Ewald_Real_set_tabulation(this, alpha)

        class(Ewald_Real), intent(inout) :: this
        real(DP), intent(in) :: alpha

        integer :: i_distance
        real(DP) :: distance_i
        
        ! cut
        do i_distance = this%i_min_distance, this%i_range_cut
            distance_i = real(i_distance, DP)*this%delta
            this%tabulation(i_distance, 1) = ewald_real_B(alpha, distance_i)
            this%tabulation(i_distance, 2) = ewald_real_C(alpha, distance_i)
        end do

        ! shift
        this%tabulation(:, 1) = this%tabulation(:, 1) - this%tabulation(this%i_range_cut, 1)
        this%tabulation(:, 2) = this%tabulation(:, 2) - this%tabulation(this%i_range_cut, 2)

    end subroutine Ewald_Real_set_tabulation

    subroutine Ewald_Real_destroy(this)
        class(Ewald_Real), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)

    end subroutine Ewald_Real_destroy

end module class_ewald_real
