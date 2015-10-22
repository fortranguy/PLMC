module types_writers_wrapper

use data_constants, only: num_components
use class_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer
use class_inter_energes_writer, only: Abstract_Inter_Energy_Writer
use class_changes_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type :: Component_Writers_Wrapper
        class(Abstract_Component_Coordinates_Writer), allocatable :: coordinates
        class(Abstract_Changes_Success_Writer), allocatable :: changes
    end type Component_Writers_Wrapper

    type, public :: Writers_Wrapper
        type(Component_Writers_Wrapper) :: components(num_components)
        class(Abstract_Inter_Energy_Writer), allocatable :: short_inter, long_inter
    end type Writers_Wrapper

end module types_writers_wrapper
