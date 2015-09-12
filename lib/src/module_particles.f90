module module_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_particles_number, only: Abstract_Particles_Number
use class_diameters, only: Abstract_Diameters
use class_moments_norm, only: Abstract_Moments_Norm
use class_positions, only: Abstract_Positions
use class_orientations, only: Abstract_Orientations

implicit none

private
public Generic_Particles_construct, Generic_Particles_destroy

    type, public :: Generic_Particle
        integer :: i_particle
        real(DP) :: diameter, moment_norm
        real(DP), dimension(num_dimensions) :: position, orientation
    end type Generic_Particle

    type, public :: Generic_Particles
        class(Abstract_Particles_Number), pointer :: particles_number
        class(Abstract_Diameters), pointer :: diameters
        class(Abstract_Moments_Norm), pointer :: moments_norm
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Orientations), pointer :: orientations
    end type Generic_Particles

contains

    subroutine Generic_Particles_construct(particles, particles_number, &
                                           diameters, moments_norm, &
                                           positions, orientations)
        type(Generic_Particles), intent(out) :: particles
        class(Abstract_Particles_Number), target, intent(in) :: particles_number
        class(Abstract_Diameters), target, intent(in) :: diameters
        class(Abstract_Moments_Norm), target, intent(in) :: moments_norm
        class(Abstract_Positions), target, intent(in) :: positions
        class(Abstract_Orientations), target, intent(in) :: orientations

        particles%particles_number => particles_number
        particles%diameters => diameters
        particles%moments_norm => moments_norm
        particles%positions => positions
        particles%orientations => orientations
    end subroutine Generic_Particles_construct

    subroutine Generic_Particles_destroy(particles)
        type(Generic_Particles), intent(inout) :: particles

        particles%orientations => null()
        particles%positions => null()
        particles%moments_norm => null()
        particles%diameters => null()
        particles%particles_number => null()
    end subroutine Generic_Particles_destroy

end module module_particles
