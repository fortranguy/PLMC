module types_particles_wrapper

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

    type, public :: Mixture_Wrapper
        type(Particles_Wrapper) :: components(2)
        class(Abstract_Particles_Diameter), allocatable :: inter_diameter
    end type Mixture_Wrapper

end module types_particles_wrapper
