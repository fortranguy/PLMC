module module_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_chemical_potential, only: Abstract_Particles_Chemical_Potential
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use class_particles_total_moment, only: Abstract_Particles_Total_Moment

implicit none

private

    type, public :: Particles_Wrapper_Parameters
        logical :: exist
        logical :: are_dipolar
        logical :: can_exchange
    end type Particles_Wrapper_Parameters

    interface assignment(=)
        module procedure :: Particles_Wrapper_Parameters_assignment
    end interface

    type, public :: Particles_Wrapper
        class(Abstract_Particles_Number), allocatable :: number
        class(Abstract_Particles_Diameter), allocatable :: diameter
        class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
        class(Abstract_Particles_Positions), allocatable :: positions
        class(Abstract_Particles_Orientations), allocatable :: orientations
        class(Abstract_Particles_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Particles_Dipolar_Moments), allocatable :: dipolar_moments
        class(Abstract_Particles_Total_Moment), allocatable :: total_moment
    end type Particles_Wrapper

contains

    pure subroutine Particles_Wrapper_Parameters_assignment(particles_parameters_target, &
        particles_parameters_value)
        type(Particles_Wrapper_Parameters), intent(out) :: particles_parameters_target
        type(Particles_Wrapper_Parameters), intent(in) :: particles_parameters_value

        particles_parameters_target%exist = particles_parameters_value%exist
        particles_parameters_target%are_dipolar = particles_parameters_value%are_dipolar
        particles_parameters_target%can_exchange = particles_parameters_value%can_exchange
    end subroutine Particles_Wrapper_Parameters_assignment

end module module_particles
