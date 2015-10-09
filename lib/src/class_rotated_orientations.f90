module class_rotated_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_increase_factor
use class_particles_orientations, only: Abstract_Particles_Orientations
use procedures_random, only: markov_orientation
use module_change_tuning, only: Concrete_Tuning_Parameters, &
    set_increase_factor
use class_changed_coordinates, only: Abstract_Changed_Coordinates

implicit none

private

    type, extends(Abstract_Changed_Coordinates), abstract, public :: Abstract_Rotated_Orientations
    private
        class(Abstract_Particles_Orientations), pointer :: orientations
        real(DP) :: delta
        type(Concrete_Tuning_Parameters) :: tuning_parameters
        real(DP) :: current_increase_factor
        logical :: max_factor_reached
    contains
        procedure :: construct => Abstract_Rotated_Orientations_construct
        procedure :: destroy => Abstract_Rotated_Orientations_destroy
        procedure :: increase_delta => Abstract_Rotated_Orientations_increase_delta
        procedure :: decrease_delta => Abstract_Rotated_Orientations_decrease_delta
        procedure :: get => Abstract_Rotated_Orientations_get
    end type Abstract_Rotated_Orientations

    type, extends(Abstract_Rotated_Orientations), public :: &
        Null_Rotated_Orientations
    contains
        procedure :: construct => Null_Rotated_Orientations_construct
        procedure :: destroy => Null_Rotated_Orientations_destroy
        procedure :: increase_delta => Null_Rotated_Orientations_increase_delta
        procedure :: decrease_delta => Null_Rotated_Orientations_decrease_delta
        procedure :: get => Null_Rotated_Orientations_get
    end type Null_Rotated_Orientations

    type, extends(Abstract_Rotated_Orientations), public :: &
        Concrete_Rotated_Orientations

    end type Concrete_Rotated_Orientations

contains

!implementation Abstract_Rotated_Orientations

    subroutine Abstract_Rotated_Orientations_construct(this, orientations, delta, &
        tuning_parameters)
        class(Abstract_Rotated_Orientations), intent(out) :: this
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations
        real(DP), intent(in) :: delta
        type(Concrete_Tuning_Parameters), intent(in) :: tuning_parameters

        this%orientations => orientations
        call check_positive("Abstract_Rotated_Orientations", "delta", delta)
        this%delta = delta
        call check_increase_factor("Abstract_Rotated_Orientations", "increase_factor", &
            tuning_parameters%increase_factor)
        this%tuning_parameters%increase_factor = tuning_parameters%increase_factor
        this%current_increase_factor = this%tuning_parameters%increase_factor
        call check_increase_factor("Abstract_Rotated_Orientations", "increase_factor_max", &
            tuning_parameters%increase_factor_max)
        this%tuning_parameters%increase_factor_max = tuning_parameters%increase_factor_max
    end subroutine Abstract_Rotated_Orientations_construct

    subroutine Abstract_Rotated_Orientations_destroy(this)
        class(Abstract_Rotated_Orientations), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_Rotated_Orientations_destroy

    subroutine Abstract_Rotated_Orientations_increase_delta(this)
        class(Abstract_Rotated_Orientations), intent(inout) :: this

        if (this%max_factor_reached) return
        call set_increase_factor("Abstract_Rotated_Orientations", this%current_increase_factor, &
            this%tuning_parameters, this%max_factor_reached)
        this%delta = this%current_increase_factor * this%delta
    end subroutine Abstract_Rotated_Orientations_increase_delta

    subroutine Abstract_Rotated_Orientations_decrease_delta(this)
        class(Abstract_Rotated_Orientations), intent(inout) :: this

        this%delta = this%delta / this%current_increase_factor
    end subroutine Abstract_Rotated_Orientations_decrease_delta

    function Abstract_Rotated_Orientations_get(this, i_particle) result(rotated_orientation)
        real(DP) :: rotated_orientation(num_dimensions)
        class(Abstract_Rotated_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle

        rotated_orientation = this%orientations%get(i_particle)
        call markov_orientation(rotated_orientation, this%delta)
    end function Abstract_Rotated_Orientations_get

!end implementation Abstract_Rotated_Orientations

!implementation Null_Rotated_Orientations

    subroutine Null_Rotated_Orientations_construct(this, orientations, delta, tuning_parameters)
        class(Null_Rotated_Orientations), intent(out) :: this
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations
        real(DP), intent(in) :: delta
        type(Concrete_Tuning_Parameters), intent(in) :: tuning_parameters
    end subroutine Null_Rotated_Orientations_construct

    subroutine Null_Rotated_Orientations_destroy(this)
        class(Null_Rotated_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Orientations_destroy

    subroutine Null_Rotated_Orientations_increase_delta(this)
        class(Null_Rotated_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Orientations_increase_delta

    subroutine Null_Rotated_Orientations_decrease_delta(this)
        class(Null_Rotated_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Orientations_decrease_delta

    function Null_Rotated_Orientations_get(this, i_particle) result(rotated_orientation)
        class(Null_Rotated_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)
        rotated_orientation = 0._DP
    end function Null_Rotated_Orientations_get

!end implementation Null_Rotated_Orientations

end module class_rotated_orientations
