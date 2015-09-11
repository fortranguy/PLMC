module class_small_rotation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use procedures_checks, only: check_positive
use class_orientations, only: Abstract_Orientations
use procedures_orientation, only: markov_orientation

implicit none

    type, abstract, public :: Abstract_Small_Rotation
    private
        class(Abstract_Orientations), pointer :: orientations
        real(DP) :: delta
    contains
        procedure :: construct => Abstract_Small_Rotation_construct
        procedure :: destroy => Abstract_Small_Rotation_destroy
        procedure :: set => Abstract_Small_Rotation_set
        procedure :: get => Abstract_Small_Rotation_get
    end type Abstract_Small_Rotation

    type, extends(Abstract_Small_Rotation), public :: Null_Small_Rotation
    contains
        procedure :: construct => Null_Small_Rotation_construct
        procedure :: destroy => Null_Small_Rotation_destroy
        procedure :: set => Null_Small_Rotation_set
        procedure :: get => Null_Small_Rotation_get
    end type Null_Small_Rotation

    type, extends(Abstract_Small_Rotation), public :: Concrete_Small_Rotation

    end type Concrete_Small_Rotation

contains

!implementation Abstract_Small_Rotation

    subroutine Abstract_Small_Rotation_construct(this, orientations)
        class(Abstract_Small_Rotation), intent(out) :: this
        class(Abstract_Orientations), target, intent(in) :: orientations

        this%orientations => orientations
    end subroutine Abstract_Small_Rotation_construct

    subroutine Abstract_Small_Rotation_destroy(this)
        class(Abstract_Small_Rotation), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_Small_Rotation_destroy

    subroutine Abstract_Small_Rotation_set(this, delta)
        class(Abstract_Small_Rotation), intent(inout) :: this
        real(DP), intent(in) :: delta

        call check_positive("Abstract_Small_Rotation", "delta", delta)
        this%delta = delta
    end subroutine Abstract_Small_Rotation_set

    function Abstract_Small_Rotation_get(this, i_particle) result(rotated_orientation)
        class(Abstract_Small_Rotation), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)

        rotated_orientation = this%orientations%get(i_particle)
        call markov_orientation(rotated_orientation, this%delta)
    end function Abstract_Small_Rotation_get

!end implementation Abstract_Small_Rotation

!implementation Null_Small_Rotation

    subroutine Null_Small_Rotation_construct(this, orientations)
        class(Null_Small_Rotation), intent(out) :: this
        class(Abstract_Orientations), target, intent(in) :: orientations
    end subroutine Null_Small_Rotation_construct

    subroutine Null_Small_Rotation_destroy(this)
        class(Null_Small_Rotation), intent(inout) :: this
    end subroutine Null_Small_Rotation_destroy

    subroutine Null_Small_Rotation_set(this, delta)
        class(Null_Small_Rotation), intent(inout) :: this
        real(DP), intent(in) :: delta
    end subroutine Null_Small_Rotation_set

    function Null_Small_Rotation_get(this, i_particle) result(rotated_orientation)
        class(Null_Small_Rotation), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: rotated_orientation(num_dimensions)
    end function Null_Small_Rotation_get

!end implementation Null_Small_Rotation

end module class_small_rotation
