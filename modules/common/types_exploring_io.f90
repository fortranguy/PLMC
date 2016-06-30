module types_exploring_io

use json_module, only: json_core, json_file, json_value
use types_exploring_readers_wrapper, only: Exploring_Readers_Wrapper
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper

implicit none

private

    type, public :: Exploring_IO_Wrapper
        type(json_core) :: json
        type(json_file) :: generating_data, exploring_data
        type(json_value), pointer :: report_data => null()
        type(Exploring_Readers_Wrapper) :: readers
        type(Exploring_Writers_Wrapper) :: writers
    end type Exploring_IO_Wrapper

end module types_exploring_io
