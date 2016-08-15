module types_component_coordinates_writer_wrapper

use classes_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer

implicit none

private

    type, public :: Component_Coordinates_Writer_Wrapper
        class(Abstract_Component_Coordinates_Writer), allocatable :: writer
    end type Component_Coordinates_Writer_Wrapper

end module types_component_coordinates_writer_wrapper
