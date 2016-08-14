module procedures_generating_readers_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_complete_coordinates_reader_factory, only: complete_coordinates_reader_create => &
    create, complete_coordinates_reader_destroy => destroy
use classes_complete_coordinates_reader, only: Abstract_Complete_Coordinates_Reader
use types_generating_readers_wrapper, only: Generating_Readers_Wrapper

implicit none

private
public :: create, destroy, set

contains

    subroutine create(readers, environment, components)
        type(Generating_Readers_Wrapper), intent(out) :: readers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)

        call complete_coordinates_reader_create(readers%complete_coordinates, environment%&
            periodic_box, components)
    end subroutine create

    subroutine destroy(readers)
        type(Generating_Readers_Wrapper), intent(inout) :: readers

        call complete_coordinates_reader_destroy(readers%complete_coordinates)
    end subroutine destroy

    subroutine set(complete_coordinates, generating_data, prefix)
        class(Abstract_Complete_Coordinates_Reader), intent(in) :: complete_coordinates
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: coordinates_filename, data_field
        logical :: data_found

        data_field = prefix//"initial coordinates"
        call generating_data%get(data_field, coordinates_filename, data_found)
        call check_data_found(data_field, data_found)
        call complete_coordinates%read(coordinates_filename)
    end subroutine set

end module procedures_generating_readers_factory
