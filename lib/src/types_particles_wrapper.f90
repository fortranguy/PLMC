module types_particles_wrapper

use data_constants, only: num_components
use class_particles_number, only: Abstract_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_component_positions, only: Abstract_Component_Positions
use class_component_orientations, only: Abstract_Component_Orientations
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use class_particles_total_moment, only: Abstract_Particles_Total_Moment

implicit none

private

    type, public :: Particles_Wrapper
        class(Abstract_Particles_Number), allocatable :: number
        class(Abstract_Particles_Diameter), allocatable :: diameter
        class(Abstract_Particles_Diameter), allocatable :: wall_diameter
        class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
        class(Abstract_Component_Positions), allocatable :: positions
        class(Abstract_Component_Orientations), allocatable :: orientations
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Particles_Dipolar_Moments), allocatable :: dipolar_moments
        class(Abstract_Particles_Total_Moment), allocatable :: total_moment
    end type Particles_Wrapper

    type, public :: Mixture_Wrapper
        type(Particles_Wrapper) :: components(num_components)
        class(Abstract_Particles_Diameter), allocatable :: inter_diameter
    end type Mixture_Wrapper

end module types_particles_wrapper
