module types_component_wrapper

use class_component_number, only: Abstract_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_component_average_number, only: Abstract_Component_Average_Number

implicit none

private

    type, public :: Component_Wrapper
        class(Abstract_Component_Number), allocatable :: number
        class(Abstract_Component_Coordinates), allocatable :: positions, orientations
        class(Abstract_Component_Dipolar_Moments), allocatable :: dipolar_moments
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Component_Average_Number), allocatable :: average_number
    end type Component_Wrapper

end module types_component_wrapper
