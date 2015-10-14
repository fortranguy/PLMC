module types_component_wrapper

use data_constants, only: num_components
use class_component_number, only: Abstract_Component_Number
use class_component_diameter, only: Abstract_Component_Diameter
use class_component_moment_norm, only: Abstract_Component_Moment_Norm
use class_component_positions, only: Abstract_Component_Positions
use class_component_orientations, only: Abstract_Component_Orientations
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_component_total_moment, only: Abstract_Component_Total_Moment

implicit none

private

    type, public :: Component_Wrapper
        class(Abstract_Component_Number), allocatable :: number
        class(Abstract_Component_Diameter), allocatable :: diameter
        class(Abstract_Component_Diameter), allocatable :: wall_diameter
        class(Abstract_Component_Moment_Norm), allocatable :: moment_norm
        class(Abstract_Component_Positions), allocatable :: positions
        class(Abstract_Component_Orientations), allocatable :: orientations
        class(Abstract_Component_Chemical_Potential), allocatable :: chemical_potential
        class(Abstract_Component_Dipolar_Moments), allocatable :: dipolar_moments
        class(Abstract_Component_Total_Moment), allocatable :: total_moment
    end type Component_Wrapper

    type, public :: Mixture_Wrapper
        type(Component_Wrapper) :: components(num_components)
        class(Abstract_Component_Diameter), allocatable :: inter_diameter
    end type Mixture_Wrapper

end module types_component_wrapper
