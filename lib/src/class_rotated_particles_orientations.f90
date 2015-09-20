module class_rotated_particles_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_positive, check_adaptation_factor
use class_particles_orientations, only: Abstract_Particles_Orientations
use procedures_orientation, only: markov_orientation

implicit none

private

    type, abstract, public :: Abstract_Rotated_Particles_Orientations
    private
        class(Abstract_Particles_Orientations), pointer :: orientations
        real(DP) :: delta
        real(DP) :: adaptation_factor
    contains
        procedure :: construct => Abstract_Rotated_Particles_Orientations_construct
        procedure :: destroy => Abstract_Rotated_Particles_Orientations_destroy
        procedure :: increase => Abstract_Rotated_Particles_Orientations_increase
        procedure :: decrease => Abstract_Rotated_Particles_Orientations_decrease
        procedure :: get => Abstract_Rotated_Particles_Orientations_get
    end type Abstract_Rotated_Particles_Orientations

    type, extends(Abstract_Rotated_Particles_Orientations), public :: &
        Null_Rotated_Particles_Orientations
    contains
        procedure :: construct => Null_Rotated_Particles_Orientations_construct
        procedure :: destroy => Null_Rotated_Particles_Orientations_destroy
        procedure :: increase => Null_Rotated_Particles_Orientations_increase
        procedure :: decrease => Null_Rotated_Particles_Orientations_decrease
        procedure :: get => Null_Rotated_Particles_Orientations_get
    end type Null_Rotated_Particles_Orientations

    type, extends(Abstract_Rotated_Particles_Orientations), public :: &
        Concrete_Rotated_Particles_Orientations

    end type Concrete_Rotated_Particles_Orientations

contains

!implementation Abstract_Rotated_Particles_Orientations

    subroutine Abstract_Rotated_Particles_Orientations_construct(this, orientations, &
        delta, adaptation_factor)
        class(Abstract_Rotated_Particles_Orientations), intent(out) :: this
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations
        real(DP), intent(in) :: delta, adaptation_factor

        this%orientations => orientations
        call check_positive("Abstract_Rotated_Particles_Orientations", "delta", delta)
        this%delta = delta
        call check_adaptation_factor("Abstract_Rotated_Particles_Orientations", &
            adaptation_factor)
        this%adaptation_factor = adaptation_factor
    end subroutine Abstract_Rotated_Particles_Orientations_construct

    subroutine Abstract_Rotated_Particles_Orientations_destroy(this)
        class(Abstract_Rotated_Particles_Orientations), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_Rotated_Particles_Orientations_destroy

    subroutine Abstract_Rotated_Particles_Orientations_increase(this)
        class(Abstract_Rotated_Particles_Orientations), intent(inout) :: this

        this%delta = this%adaptation_factor * this%delta
    end subroutine Abstract_Rotated_Particles_Orientations_increase

    subroutine Abstract_Rotated_Particles_Orientations_decrease(this)
        class(Abstract_Rotated_Particles_Orientations), intent(inout) :: this

        this%delta = this%delta / this%adaptation_factor
    end subroutine Abstract_Rotated_Particles_Orientations_decrease

    function Abstract_Rotated_Particles_Orientations_get(this, i_particle) &
        result(rotated_orientation)
        class(Abstract_Rotated_Particles_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)

        rotated_orientation = this%orientations%get(i_particle)
        call markov_orientation(rotated_orientation, this%delta)
    end function Abstract_Rotated_Particles_Orientations_get

!end implementation Abstract_Rotated_Particles_Orientations

!implementation Null_Rotated_Particles_Orientations

    subroutine Null_Rotated_Particles_Orientations_construct(this, orientations, delta, &
        adaptation_factor)
        class(Null_Rotated_Particles_Orientations), intent(out) :: this
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations
        real(DP), intent(in) :: delta, adaptation_factor
    end subroutine Null_Rotated_Particles_Orientations_construct

    subroutine Null_Rotated_Particles_Orientations_destroy(this)
        class(Null_Rotated_Particles_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Particles_Orientations_destroy

    subroutine Null_Rotated_Particles_Orientations_increase(this)
        class(Null_Rotated_Particles_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Particles_Orientations_increase

    subroutine Null_Rotated_Particles_Orientations_decrease(this)
        class(Null_Rotated_Particles_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Particles_Orientations_decrease

    function Null_Rotated_Particles_Orientations_get(this, i_particle) &
        result(rotated_orientation)
        class(Null_Rotated_Particles_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)
        rotated_orientation = 0._DP
    end function Null_Rotated_Particles_Orientations_get

!end implementation Null_Rotated_Particles_Orientations

end module class_rotated_particles_orientations
