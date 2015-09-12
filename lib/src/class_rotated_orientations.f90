module class_rotated_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_positive
use class_orientations, only: Abstract_Orientations
use procedures_orientation, only: markov_orientation

implicit none

    type, abstract, public :: Abstract_Rotated_Orientations
    private
        class(Abstract_Orientations), pointer :: orientations
        real(DP) :: delta
    contains
        procedure :: construct => Abstract_Rotated_Orientations_construct
        procedure :: destroy => Abstract_Rotated_Orientations_destroy
        procedure :: set => Abstract_Rotated_Orientations_set
        procedure :: get => Abstract_Rotated_Orientations_get
    end type Abstract_Rotated_Orientations

    type, extends(Abstract_Rotated_Orientations), public :: Null_Rotated_Orientations
    contains
        procedure :: construct => Null_Rotated_Orientations_construct
        procedure :: destroy => Null_Rotated_Orientations_destroy
        procedure :: set => Null_Rotated_Orientations_set
        procedure :: get => Null_Rotated_Orientations_get
    end type Null_Rotated_Orientations

    type, extends(Abstract_Rotated_Orientations), public :: Concrete_Rotated_Orientations

    end type Concrete_Rotated_Orientations

contains

!implementation Abstract_Rotated_Orientations

    subroutine Abstract_Rotated_Orientations_construct(this, orientations)
        class(Abstract_Rotated_Orientations), intent(out) :: this
        class(Abstract_Orientations), target, intent(in) :: orientations

        this%orientations => orientations
    end subroutine Abstract_Rotated_Orientations_construct

    subroutine Abstract_Rotated_Orientations_destroy(this)
        class(Abstract_Rotated_Orientations), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_Rotated_Orientations_destroy

    subroutine Abstract_Rotated_Orientations_set(this, delta)
        class(Abstract_Rotated_Orientations), intent(inout) :: this
        real(DP), intent(in) :: delta

        call check_positive("Abstract_Rotated_Orientations", "delta", delta)
        this%delta = delta
    end subroutine Abstract_Rotated_Orientations_set

    function Abstract_Rotated_Orientations_get(this, i_particle) result(rotated_orientation)
        class(Abstract_Rotated_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)

        rotated_orientation = this%orientations%get(i_particle)
        call markov_orientation(rotated_orientation, this%delta)
    end function Abstract_Rotated_Orientations_get

!end implementation Abstract_Rotated_Orientations

!implementation Null_Rotated_Orientations

    subroutine Null_Rotated_Orientations_construct(this, orientations)
        class(Null_Rotated_Orientations), intent(out) :: this
        class(Abstract_Orientations), target, intent(in) :: orientations
    end subroutine Null_Rotated_Orientations_construct

    subroutine Null_Rotated_Orientations_destroy(this)
        class(Null_Rotated_Orientations), intent(inout) :: this
    end subroutine Null_Rotated_Orientations_destroy

    subroutine Null_Rotated_Orientations_set(this, delta)
        class(Null_Rotated_Orientations), intent(inout) :: this
        real(DP), intent(in) :: delta
    end subroutine Null_Rotated_Orientations_set

    function Null_Rotated_Orientations_get(this, i_particle) result(rotated_orientation)
        class(Null_Rotated_Orientations), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)
        rotated_orientation = 0._DP
    end function Null_Rotated_Orientations_get

!end implementation Null_Rotated_Orientations

end module class_rotated_orientations
