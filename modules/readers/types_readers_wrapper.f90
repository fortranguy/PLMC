module types_readers_wrapper

use classes_component_coordinates_reader, only: Abstract_Component_Coordinates_Reader

implicit none

private

    type, public :: Component_Coordinates_Reader_wrapper
        class(Abstract_Component_Coordinates_Reader), allocatable :: coordinates
    end type Component_Coordinates_Reader_wrapper

    type, public :: Readers_Wrapper
        type(Component_Coordinates_Reader_wrapper), allocatable :: components(:)
    end type Readers_Wrapper

end module types_readers_wrapper
