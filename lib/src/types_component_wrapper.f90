module types_component_wrapper

use data_constants, only: num_components
use class_component_number, only: Abstract_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_component_total_moment, only: Abstract_Component_Total_Moment

implicit none

private

    type, public :: Component_Wrapper
        class(Abstract_Component_Number), allocatable :: number
        class(Abstract_Component_Coordinates), allocatable :: positions, orientations
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Component_Dipolar_Moments), allocatable :: dipolar_moments
        class(Abstract_Component_Total_Moment), allocatable :: total_moment
    end type Component_Wrapper

    type, public :: Mixture_Wrapper_Old
        type(Component_Wrapper) :: components(num_components)
    end type Mixture_Wrapper_Old

end module types_component_wrapper
