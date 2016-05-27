module classes_component_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use classes_coordinates, only: Abstract_Coordinates
use classes_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Component_Dipolar_Moments
    private
        class(Abstract_Component_Coordinates), pointer :: orientations => null()
        real(DP) :: norm = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num => Abstract_get_num
        procedure :: get_norm => Abstract_get_norm
        procedure :: get => Abstract_get
    end type Abstract_Component_Dipolar_Moments

    type, extends(Abstract_Component_Dipolar_Moments), public :: Concrete_Component_Dipolar_Moments

    end type Concrete_Component_Dipolar_Moments

    type, extends(Abstract_Component_Dipolar_Moments), public :: Null_Component_Dipolar_Moments
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_num => Null_get_num
        procedure :: get_norm => Null_get_norm
        procedure :: get => Null_get
    end type Null_Component_Dipolar_Moments

contains

!implementation Abstract_Component_Dipolar_Moments

    subroutine Abstract_construct(this, norm, orientations)
        class(Abstract_Component_Dipolar_Moments), intent(out) :: this
        real(DP), intent(in) :: norm
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations

        call check_positive("Abstract_construct", "norm", norm)
        this%norm = norm
        this%orientations => orientations
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Component_Dipolar_Moments), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num(this) result(num_moments)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this

        num_moments = this%orientations%get_num()
    end function Abstract_get_num

    pure real(DP) function Abstract_get_norm(this) result(norm)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this

        norm = this%norm
    end function Abstract_get_norm

    pure function Abstract_get(this, i_particle) result(dipolar_moment)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)

        dipolar_moment = this%norm * this%orientations%get(i_particle)
    end function Abstract_get

!end implementation Abstract_Component_Dipolar_Moments

!implementation Null_Component_Dipolar_Moments

    subroutine Null_construct(this, norm, orientations)
        class(Null_Component_Dipolar_Moments), intent(out) :: this
        real(DP), intent(in) :: norm
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Component_Dipolar_Moments), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get_num(this) result(num_moments)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        num_moments = 0
    end function Null_get_num

    pure real(DP) function Null_get_norm(this) result(norm)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        norm = 0._DP
    end function Null_get_norm

    pure function Null_get(this, i_particle) result(dipolar_moment)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)
        dipolar_moment = 0._DP
    end function Null_get

!end implementation Null_Component_Dipolar_Moments

end module classes_component_dipolar_moments
