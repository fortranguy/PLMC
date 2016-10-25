module procedures_readers_factory

use data_input_prefixes, only: readers_prefix
use json_module, only: json_file
use procedures_checks, only: check_data_found
use procedures_command_arguments, only: create_filename_from_argument
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_complete_coordinates_reader_factory, only: complete_coordinates_reader_create => &
    create, complete_coordinates_reader_destroy => destroy
use types_readers_wrapper, only: Readers_Wrapper

implicit none

private
public :: create, destroy, set

interface set
    module procedure :: set_from_json
    module procedure :: set_from_snap
end interface

contains

    subroutine create(readers, environment, components, particle_insertion_domains)
        type(Readers_Wrapper), intent(out) :: readers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:, :)
        class(Abstract_Parallelepiped_Domain), optional, intent(in) ::particle_insertion_domains(:)

        call complete_coordinates_reader_create(readers%complete_coordinates, environment, &
            components, particle_insertion_domains)
    end subroutine create

    subroutine destroy(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call complete_coordinates_reader_destroy(readers%complete_coordinates)
    end subroutine destroy

    !> @note var_type is optional, it shouldn't be necessary.
    !> @todo send a bug report to gfortran and ifort?
    subroutine set_from_json(readers, generating_data)
        type(Readers_Wrapper), intent(inout) :: readers
        type(json_file), intent(inout) :: generating_data

        type(String_Wrapper), allocatable :: coordinates(:)
        type(Concrete_Number_to_String) :: string
        integer :: dummy_var_type, num_boxes, i_box
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = readers_prefix//"initial coordinates"
        call generating_data%info(data_field, data_found, dummy_var_type, num_boxes)
        call check_data_found(data_field, data_found)
        allocate(coordinates(num_boxes))
        do i_box = 1, size(coordinates)
            call generating_data%get(data_field//"["//string%get(i_box)//"]", coordinates(i_box)%&
                string, data_found)
            call check_data_found(data_field, data_found)
        end do
        call readers%complete_coordinates%read(coordinates)
    end subroutine set_from_json

    subroutine set_from_snap(readers, offset_i_snap)
        type(Readers_Wrapper), intent(inout) :: readers
        integer, intent(in) :: offset_i_snap

        type(String_Wrapper) :: coordinates(1)

        call create_filename_from_argument(coordinates(1)%string, offset_i_snap)
        call readers%complete_coordinates%read(coordinates)
    end subroutine set_from_snap

end module procedures_readers_factory
