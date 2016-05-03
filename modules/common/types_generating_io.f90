module types_generating_io

use json_module, only: json_file
use types_generating_readers_wrapper, only: Generating_Readers_Wrapper
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper

implicit none

private

    type, public :: Generating_IO_Wrapper
        type(json_file) :: generating_data
        type(Generating_Readers_Wrapper) :: readers
        type(Generating_Writers_Wrapper) :: writers
    end type Generating_IO_Wrapper

end module types_generating_io
