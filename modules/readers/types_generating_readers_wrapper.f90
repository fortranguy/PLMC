module types_generating_readers_wrapper

use classes_complete_coordinates_reader, only: Abstract_Complete_Coordinates_Reader

implicit none

private

    type, public :: Generating_Readers_Wrapper
        class(Abstract_Complete_Coordinates_Reader), allocatable :: complete_coordinates
    end type Generating_Readers_Wrapper

end module types_generating_readers_wrapper
