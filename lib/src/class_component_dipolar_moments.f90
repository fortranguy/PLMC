module class_component_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use class_coordinates, only: Abstract_Coordinates
use class_component_coordinates, only: Abstract_Component_Coordinates

implicit none

private

    type, extends(Abstract_Coordinates), abstract, public :: Abstract_Component_Dipolar_Moments
    private
        real(DP) :: norm
        class(Abstract_Component_Coordinates), pointer :: orientations
    contains
        procedure :: construct => Abstract_Component_Dipolar_Moments_construct
        procedure :: destroy => Abstract_Component_Dipolar_Moments_destroy
        procedure :: get_num => Abstract_Component_Dipolar_Moments_get_num
        procedure :: get => Abstract_Component_Dipolar_Moments_get
        procedure :: get_norm => Abstract_Component_Dipolar_Moments_get_norm
    end type Abstract_Component_Dipolar_Moments

    type, extends(Abstract_Component_Dipolar_Moments), public :: Concrete_Component_Dipolar_Moments

    end type Concrete_Component_Dipolar_Moments

    type, extends(Abstract_Component_Dipolar_Moments), public :: Null_Component_Dipolar_Moments
    contains
        procedure :: construct => Null_Component_Dipolar_Moments_construct
        procedure :: destroy => Null_Component_Dipolar_Moments_destroy
        procedure :: get_num => Null_Component_Dipolar_Moments_get_num
        procedure :: get => Null_Component_Dipolar_Moments_get
        procedure :: get_norm => Null_Component_Dipolar_Moments_get_norm
    end type Null_Component_Dipolar_Moments

contains

!implementation Abstract_Component_Dipolar_Moments

    subroutine Abstract_Component_Dipolar_Moments_construct(this, norm, orientations)
        class(Abstract_Component_Dipolar_Moments), intent(out) :: this
        real(DP), intent(in) :: norm
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations

        call check_positive("Abstract_Component_Dipolar_Moments_construct", "norm", norm)
        this%norm = norm
        this%orientations => orientations
    end subroutine Abstract_Component_Dipolar_Moments_construct

    subroutine Abstract_Component_Dipolar_Moments_destroy(this)
        class(Abstract_Component_Dipolar_Moments), intent(inout) :: this

        this%orientations => null()
    end subroutine Abstract_Component_Dipolar_Moments_destroy

    pure integer function Abstract_Component_Dipolar_Moments_get_num(this) result(num_moments)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this

        num_moments = this%orientations%get_num()
    end function Abstract_Component_Dipolar_Moments_get_num

    pure function Abstract_Component_Dipolar_Moments_get(this, i_particle) result(moment)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: moment(num_dimensions)

        moment = this%norm * this%orientations%get(i_particle)
    end function Abstract_Component_Dipolar_Moments_get

    pure real(DP) function Abstract_Component_Dipolar_Moments_get_norm(this) result(norm)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: this

        norm = this%norm
    end function Abstract_Component_Dipolar_Moments_get_norm

!end implementation Abstract_Component_Dipolar_Moments

!implementation Null_Component_Dipolar_Moments

    subroutine Null_Component_Dipolar_Moments_construct(this, norm, orientations)
        class(Null_Component_Dipolar_Moments), intent(out) :: this
        real(DP), intent(in) :: norm
        class(Abstract_Component_Coordinates), target, intent(in) :: orientations
    end subroutine Null_Component_Dipolar_Moments_construct

    subroutine Null_Component_Dipolar_Moments_destroy(this)
        class(Null_Component_Dipolar_Moments), intent(inout) :: this
    end subroutine Null_Component_Dipolar_Moments_destroy

    pure integer function Null_Component_Dipolar_Moments_get_num(this) result(num_moments)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        num_moments = 0
    end function Null_Component_Dipolar_Moments_get_num

    pure function Null_Component_Dipolar_Moments_get(this, i_particle) result(moment)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: moment(num_dimensions)
        moment = 0._DP
    end function Null_Component_Dipolar_Moments_get

    pure real(DP) function Null_Component_Dipolar_Moments_get_norm(this) result(norm)
        class(Null_Component_Dipolar_Moments), intent(in) :: this
        norm = 0._DP
    end function Null_Component_Dipolar_Moments_get_norm

!end implementation Null_Component_Dipolar_Moments

end module class_component_dipolar_moments
