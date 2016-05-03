module procedures_readers_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use types_component_wrapper, only: Component_Wrapper
use procedures_coordinates_reader_factory, only: coordinates_reader_create => create, &
    coordinates_reader_destroy => destroy
use types_readers_wrapper, only: Component_Coordinates_Reader_wrapper, Readers_Wrapper

implicit none

private
public :: readers_create, readers_destroy, readers_set_initial_coordinates

contains

    subroutine readers_create(readers, components)
        type(Readers_Wrapper), intent(out) :: readers
        type(Component_Wrapper), intent(in) :: components(:)

        call coordinates_reader_create(readers%components, components)
    end subroutine readers_create

    subroutine readers_destroy(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call coordinates_reader_destroy(readers%components)
    end subroutine readers_destroy


    subroutine readers_set_initial_coordinates(components_readers, input_data, prefix)
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
    end subroutine readers_set_initial_coordinates

end module procedures_readers_factory
