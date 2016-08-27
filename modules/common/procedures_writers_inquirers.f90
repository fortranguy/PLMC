module procedures_writers_inquirers

use json_module, only: json_file
use procedures_property_inquirers, only: logical_from_json

implicit none

private
public :: write_coordinates

interface write_coordinates
    module procedure :: write_coordinates_from_json
end interface write_coordinates

contains

    logical function write_coordinates_from_json(generating_data, prefix) result(write_coordinates)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        write_coordinates = logical_from_json(generating_data, prefix//"Coordinates.write")
    end function write_coordinates_from_json

end module procedures_writers_inquirers
