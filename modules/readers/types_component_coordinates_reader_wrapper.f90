module types_component_coordinates_reader_wrapper

use classes_component_coordinates_reader, only: Abstract_Component_Coordinates_Reader

implicit none

private

    type, public :: Component_Coordinates_Reader_wrapper
        class(Abstract_Component_Coordinates_Reader), allocatable :: reader
    end type Component_Coordinates_Reader_wrapper

end module types_component_coordinates_reader_wrapper
