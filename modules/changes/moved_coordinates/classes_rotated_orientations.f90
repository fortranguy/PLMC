module classes_rotated_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_increase_factor
use classes_component_coordinates, only: Abstract_Component_Coordinates
use procedures_random_number, only: markov_orientation
use module_move_tuning, only: Concrete_Move_Tuning_Parameters, set_increase_factor
use classes_moved_coordinates, only: Abstract_Moved_Coordinates

implicit none

private

    type, extends(Abstract_Moved_Coordinates), public :: Concrete_Rotated_Orientations
    private
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        real(DP) :: delta = 0._DP
        type(Concrete_Move_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor = 1._DP
        logical :: max_factor_reached = .false.
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: increase_delta => Concrete_increase_delta
        procedure :: decrease_delta => Concrete_decrease_delta
        procedure :: get => Concrete_get
    end type Concrete_Rotated_Orientations

contains

    subroutine Concrete_construct(this, orientations, delta, tuning_parameters)
        class(Concrete_Rotated_Orientations), intent(out) :: this
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations
        real(DP), intent(in) :: delta
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters

        this%orientations => orientations
        call check_positive("Concrete_Rotated_Orientations", "delta", delta)
        this%delta = delta
        call check_increase_factor("Concrete_Rotated_Orientations", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Concrete_Rotated_Orientations", "increase_factor_max", &
            tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Concrete_Rotated_Orientations), intent(inout) :: this

        this%orientations => null()
    end subroutine Concrete_destroy

    subroutine Concrete_increase_delta(this)
        class(Concrete_Rotated_Orientations), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Concrete_Rotated_Orientations", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Concrete_increase_delta

    subroutine Concrete_decrease_delta(this)
        class(Concrete_Rotated_Orientations), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Concrete_decrease_delta

    function Concrete_get(this, i_particle) result(rotated_orientation)
        real(DP) :: rotated_orientation(num_dimensions)
        class(Concrete_Rotated_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle

        rotated_orientation = this%orientations%get(i_particle)
        call markov_orientation(rotated_orientation, this%delta)
    end function Concrete_get

end module classes_rotated_orientations
