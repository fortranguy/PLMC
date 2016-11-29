module types_exploring_io

use json_module, only: json_file
use types_json_report, only: Exploring_JSON_Report
use types_readers_wrapper, only: Readers_Wrapper
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper

implicit none

private

    type, public :: Exploring_IO_Wrapper
        type(json_file) :: data
        type(Exploring_JSON_Report) :: report
        type(Readers_Wrapper) :: readers
        type(Exploring_Writers_Wrapper) :: writers
    end type Exploring_IO_Wrapper

end module types_exploring_io
