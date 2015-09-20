module class_particles_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_orientations, only: Abstract_Particles_Orientations

implicit none

private

    type, public :: Particles_Dipolar_Moments_Facade
    private
        class(Abstract_Particles_Moment_Norm), pointer :: moment_norm
        class(Abstract_Particles_Orientations), pointer :: orientations
    contains
        procedure :: construct => Particles_Dipolar_Moments_Facade_construct
        procedure :: destroy => Particles_Dipolar_Moments_Facade_destroy
        procedure :: get_num => Particles_Dipolar_Moments_Facade_get_num
        procedure :: get => Particles_Dipolar_Moments_Facade_get
    end type Particles_Dipolar_Moments_Facade

contains

    subroutine Particles_Dipolar_Moments_Facade_construct(this, moment_norm, orientations)
        class(Particles_Dipolar_Moments_Facade), intent(out) :: this
        class(Abstract_Particles_Moment_Norm), target, intent(in) :: moment_norm
        class(Abstract_Particles_Orientations), target, intent(in) :: orientations

        this%moment_norm => moment_norm
        this%orientations => orientations
    end subroutine Particles_Dipolar_Moments_Facade_construct

    subroutine Particles_Dipolar_Moments_Facade_destroy(this)
        class(Particles_Dipolar_Moments_Facade), intent(inout) :: this

        this%orientations => null()
        this%moment_norm => null()
    end subroutine Particles_Dipolar_Moments_Facade_destroy

    pure function Particles_Dipolar_Moments_Facade_get_num(this) result(num_dipolar_moments)
        class(Particles_Dipolar_Moments_Facade), intent(in) :: this
        integer :: num_dipolar_moments

        num_dipolar_moments = this%orientations%get_num()
    end function Particles_Dipolar_Moments_Facade_get_num

    pure function Particles_Dipolar_Moments_Facade_get(this, i_particle) result(dipolar_moment)
        class(Particles_Dipolar_Moments_Facade), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: dipolar_moment(num_dimensions)

        dipolar_moment = this%moment_norm%get() * this%orientations%get(i_particle)
    end function Particles_Dipolar_Moments_Facade_get

end module class_particles_dipolar_moments
