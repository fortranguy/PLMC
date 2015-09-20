module types_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade
use class_particles_total_moment, only: Particles_Total_Moment_Facade

implicit none

private

    type, public :: Concrete_Particle
        logical :: same_type = .false.
        integer :: i = 0
        real(DP) :: diameter = 0._DP
        real(DP) :: min_diameter = 0._DP
        real(DP) :: moment_norm = 0._DP
        real(DP) :: position(num_dimensions) = 0._DP
        real(DP) :: orientation(num_dimensions) = 0._DP
    end type Concrete_Particle

    type, public :: Concrete_Particles
        class(Abstract_Particles_Number), allocatable :: number
        class(Abstract_Particles_Diameter), allocatable :: diameter
        class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
        class(Abstract_Particles_Positions), allocatable :: positions
        class(Abstract_Particles_Orientations), allocatable :: orientations
        class(Particles_Dipolar_Moments_Facade), allocatable :: dipolar_moments
        class(Particles_Total_Moment_Facade), allocatable :: total_moment
    end type Concrete_Particles

end module types_particles
