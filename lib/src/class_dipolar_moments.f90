module class_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_moment_norm, only: Abstract_Moment_Norm
use class_orientations, only: Abstract_Orientations

implicit none

private

    type, public :: Dipolar_Moments_Facade
    private
        class(Abstract_Moment_Norm), pointer :: moment_norm
        class(Abstract_Orientations), pointer :: orientations
    contains
        procedure :: construct => Dipolar_Moments_Facade_construct
        procedure :: destroy => Dipolar_Moments_Facade_destroy
        procedure :: get_num => Dipolar_Moments_Facade_get_num
        procedure :: get => Dipolar_Moments_Facade_get
    end type Dipolar_Moments_Facade

contains

    subroutine Dipolar_Moments_Facade_construct(this, moment_norm, orientations)
        class(Dipolar_Moments_Facade), intent(out) :: this
        class(Abstract_Moment_Norm), target, intent(in) :: moment_norm
        class(Abstract_Orientations), target, intent(in) :: orientations

        this%moment_norm => moment_norm
        this%orientations => orientations
    end subroutine Dipolar_Moments_Facade_construct

    subroutine Dipolar_Moments_Facade_destroy(this)
        class(Dipolar_Moments_Facade), intent(inout) :: this

        this%orientations => null()
        this%moment_norm => null()
    end subroutine Dipolar_Moments_Facade_destroy

    pure function Dipolar_Moments_Facade_get_num(this) result(num_dipolar_moments)
        class(Dipolar_Moments_Facade), intent(in) :: this
        integer :: num_dipolar_moments

        num_dipolar_moments = this%orientations%get_num()
    end function Dipolar_Moments_Facade_get_num

    pure function Dipolar_Moments_Facade_get(this, i_particle) result(dipolar_moment)
        class(Dipolar_Moments_Facade), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)

        dipolar_moment = this%moment_norm%get() * this%orientations%get(i_particle)
    end function Dipolar_Moments_Facade_get

end module class_dipolar_moments
