module types_component_wrapper

use classes_component_number, only: Abstract_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_component_average_number, only: Abstract_Component_Average_Number

implicit none

private

    type, public :: Component_Wrapper
        class(Abstract_Component_Number), allocatable :: number
        class(Abstract_Component_Coordinates), allocatable :: positions, orientations
        class(Abstract_Component_Dipole_Moments), allocatable :: dipole_moments
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Component_Average_Number), allocatable :: average_number
    end type Component_Wrapper

end module types_component_wrapper
