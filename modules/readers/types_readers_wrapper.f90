module types_readers_wrapper

use classes_box_size_reader, only: Abstract_Box_Size_Reader
use classes_component_coordinates_reader, only: Abstract_Coordinates_Reader

implicit none

private

    type, public :: Component_Readers_wrapper
        class(Abstract_Coordinates_Reader), allocatable :: coordinates
    end type Component_Readers_wrapper

    type, public :: Readers_Wrapper
        class(Abstract_Box_Size_Reader), allocatable :: box_size
        type(Component_Readers_wrapper), allocatable :: components(:)
    end type Readers_Wrapper

end module types_readers_wrapper
