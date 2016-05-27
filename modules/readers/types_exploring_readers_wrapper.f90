module types_exploring_readers_wrapper

use classes_box_size_reader, only: Abstract_Box_Size_Reader
use classes_snap_number, only: Abstract_Snap_Number
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper

implicit none

private

    type, public :: Exploring_Readers_Wrapper
        class(Abstract_Box_Size_Reader), allocatable :: box
        class(Abstract_Snap_Number), allocatable :: snap_number
        type(Component_Coordinates_Reader_wrapper), allocatable :: components(:)
    end type Exploring_Readers_Wrapper

end module types_exploring_readers_wrapper
