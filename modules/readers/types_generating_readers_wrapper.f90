module types_generating_readers_wrapper

use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper

implicit none

private

    type, public :: Generating_Readers_Wrapper
        type(Component_Coordinates_Reader_wrapper), allocatable :: components(:)
    end type Generating_Readers_Wrapper

end module types_generating_readers_wrapper
