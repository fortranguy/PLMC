module class_particles_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_orientations, only: Abstract_Particles_Orientations

implicit none

private

    type, abstract, public :: Abstract_Particles_Dipolar_Moments
    private
        class(Abstract_Particles_Moment_Norm), pointer :: moment_norm
        class(Abstract_Particles_Orientations), pointer :: orientations
    contains
        procedure :: construct => Abstract_Particles_Dipolar_Moments_construct
        procedure :: destroy => Abstract_Particles_Dipolar_Moments_destroy
        procedure :: get_num => Abstract_Particles_Dipolar_Moments_get_num
        procedure :: get => Abstract_Particles_Dipolar_Moments_get
    end type Abstract_Particles_Dipolar_Moments

    type, extends(Abstract_Particles_Dipolar_Moments), public :: Concrete_Particles_Dipolar_Moments

    end type Concrete_Particles_Dipolar_Moments

    type, extends(Abstract_Particles_Dipolar_Moments), public :: Null_Particles_Dipolar_Moments
    contains
        procedure :: construct => Null_Particles_Dipolar_Moments_construct
        procedure :: destroy => Null_Particles_Dipolar_Moments_destroy
        procedure :: get_num => Null_Particles_Dipolar_Moments_get_num
        procedure :: get => Null_Particles_Dipolar_Moments_get
    end type Null_Particles_Dipolar_Moments

contains

!implementation Abstract_Particles_Dipolar_Moments

    subroutine Abstract_Particles_Dipolar_Moments_construct(this, moment_norm, orientations)
        class(Abstract_Particles_Dipolar_Moments), intent(out) :: this
        class(Abstract_Particles_Moment_Norm), target, intent(in) :: moment_norm
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations

        this%moment_norm => moment_norm
        this%orientations => orientations
    end subroutine Abstract_Particles_Dipolar_Moments_construct

    subroutine Abstract_Particles_Dipolar_Moments_destroy(this)
        class(Abstract_Particles_Dipolar_Moments), intent(inout) :: this

        this%orientations => null()
        this%moment_norm => null()
    end subroutine Abstract_Particles_Dipolar_Moments_destroy

    pure function Abstract_Particles_Dipolar_Moments_get_num(this) result(num_dipolar_moments)
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: this
        integer :: num_dipolar_moments

        num_dipolar_moments = this%orientations%get_num()
    end function Abstract_Particles_Dipolar_Moments_get_num

    pure function Abstract_Particles_Dipolar_Moments_get(this, i_particle) result(dipolar_moment)
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)

        dipolar_moment = this%moment_norm%get() * this%orientations%get(i_particle)
    end function Abstract_Particles_Dipolar_Moments_get

!end implementation Abstract_Particles_Dipolar_Moments

!implementation Null_Particles_Dipolar_Moments

    subroutine Null_Particles_Dipolar_Moments_construct(this, moment_norm, orientations)
        class(Null_Particles_Dipolar_Moments), intent(out) :: this
        class(Abstract_Particles_Moment_Norm), target, intent(in) :: moment_norm
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations
    end subroutine Null_Particles_Dipolar_Moments_construct

    subroutine Null_Particles_Dipolar_Moments_destroy(this)
        class(Null_Particles_Dipolar_Moments), intent(inout) :: this
    end subroutine Null_Particles_Dipolar_Moments_destroy

    pure function Null_Particles_Dipolar_Moments_get_num(this) result(num_dipolar_moments)
        class(Null_Particles_Dipolar_Moments), intent(in) :: this
        integer :: num_dipolar_moments
        num_dipolar_moments = 0
    end function Null_Particles_Dipolar_Moments_get_num

    pure function Null_Particles_Dipolar_Moments_get(this, i_particle) result(dipolar_moment)
        class(Null_Particles_Dipolar_Moments), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)
        dipolar_moment = 0._DP
    end function Null_Particles_Dipolar_Moments_get

!end implementation Null_Particles_Dipolar_Moments

end module class_particles_dipolar_moments
