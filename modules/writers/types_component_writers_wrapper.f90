module types_component_writers_wrapper

use classes_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer
use classes_changes_success_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type, public :: Component_Writers_Wrapper
        class(Abstract_Component_Coordinates_Writer), allocatable :: coordinates
        class(Abstract_Changes_Success_Writer), allocatable :: changes
    end type Component_Writers_Wrapper

end module types_component_writers_wrapper
