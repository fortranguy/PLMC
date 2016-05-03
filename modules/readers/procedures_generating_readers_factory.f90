module procedures_generating_readers_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use types_component_wrapper, only: Component_Wrapper
use procedures_coordinates_reader_factory, only: coordinates_reader_create => create, &
    coordinates_reader_destroy => destroy
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper
use types_generating_readers_wrapper, only: Generating_Readers_Wrapper

implicit none

private
public :: generating_readers_create, generating_readers_destroy, generating_readers_set

contains

    subroutine generating_readers_create(readers, components)
        type(Generating_Readers_Wrapper), intent(out) :: readers
        type(Component_Wrapper), intent(in) :: components(:)

        call coordinates_reader_create(readers%components, components)
    end subroutine generating_readers_create

    subroutine generating_readers_destroy(readers)
        type(Generating_Readers_Wrapper), intent(inout) :: readers

        call coordinates_reader_destroy(readers%components)
    end subroutine generating_readers_destroy

    subroutine generating_readers_set(components_readers, input_data, prefix)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_readers(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field, filename
        logical :: data_found

        do i_component = 1, size(components_readers)
            data_field = prefix//"Component "//string%get(i_component)//"."//"initial coordinates"
            call input_data%get(data_field, filename, data_found)
            call check_data_found(data_field, data_found)
            call components_readers(i_component)%coordinates%read(filename)
        end do
    end subroutine generating_readers_set

end module procedures_generating_readers_factory
