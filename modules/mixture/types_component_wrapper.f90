module types_component_wrapper

use classes_num_particles, only: Abstract_Num_Particles
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments

implicit none

private

    type, public :: Component_Wrapper
        class(Abstract_Num_Particles), allocatable :: num_particles
        class(Abstract_Component_Coordinates), allocatable :: positions, orientations
        class(Abstract_Component_Dipole_Moments), allocatable :: dipole_moments
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
    end type Component_Wrapper

end module types_component_wrapper
