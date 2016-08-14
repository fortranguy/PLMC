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
public :: create, destroy, set

contains

    subroutine create(readers, components)
        type(Generating_Readers_Wrapper), intent(out) :: readers
        type(Component_Wrapper), intent(in) :: components(:)

        call coordinates_reader_create(readers%components, components)
    end subroutine create

    subroutine destroy(readers)
        type(Generating_Readers_Wrapper), intent(inout) :: readers

        call coordinates_reader_destroy(readers%components)
    end subroutine destroy

    subroutine set(components, generating_data, prefix)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field, filename
        logical :: data_found

        do i_component = 1, size(components)
            data_field = prefix//"Component "//string%get(i_component)//"."//"initial coordinates"
            call generating_data%get(data_field, filename, data_found)
            call check_data_found(data_field, data_found)
            !call components(i_component)%coordinates%read(filename)
        end do
    end subroutine set

end module procedures_generating_readers_factory
