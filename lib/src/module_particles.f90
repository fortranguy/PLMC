module module_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_number, only: Abstract_Number
use class_diameter, only: Abstract_Diameter
use class_moment_norm, only: Abstract_Moment_Norm
use class_positions, only: Abstract_Positions
use class_orientations, only: Abstract_Orientations

implicit none

private
public Concrete_Particles_construct, Concrete_Particles_destroy

    type, public :: Concrete_Particle
        logical :: same_type
        integer :: i
        real(DP) :: diameter, min_diameter, moment_norm
        real(DP), dimension(num_dimensions) :: position, orientation
    end type Concrete_Particle

    type, public :: Concrete_Particles
        class(Abstract_Number), pointer :: number
        class(Abstract_Diameter), pointer :: diameter
        class(Abstract_Moment_Norm), pointer :: moment_norm
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Orientations), pointer :: orientations
    end type Concrete_Particles

contains

    subroutine Concrete_Particles_construct(particles, number, &
                                           diameter, moment_norm, &
                                           positions, orientations)
        type(Concrete_Particles), intent(out) :: particles
        class(Abstract_Number), target, intent(in) :: number
        class(Abstract_Diameter), target, intent(in) :: diameter
        class(Abstract_Moment_Norm), target, intent(in) :: moment_norm
        class(Abstract_Positions), target, intent(in) :: positions
        class(Abstract_Orientations), target, intent(in) :: orientations

        particles%number => number
        particles%diameter => diameter
        particles%moment_norm => moment_norm
        particles%positions => positions
        particles%orientations => orientations
    end subroutine Concrete_Particles_construct

    subroutine Concrete_Particles_destroy(particles)
        type(Concrete_Particles), intent(inout) :: particles

        particles%orientations => null()
        particles%positions => null()
        particles%moment_norm => null()
        particles%diameter => null()
        particles%number => null()
    end subroutine Concrete_Particles_destroy

end module module_particles
