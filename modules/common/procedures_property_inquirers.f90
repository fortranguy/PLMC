module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found


implicit none

private
public :: logical_from_json

contains

    logical function logical_from_json(input_data, statement)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: statement

        logical :: data_found

        call input_data%get(statement, logical_from_json, data_found)
        call check_data_found(statement, data_found)
    end function logical_from_json

end module procedures_property_inquirers
