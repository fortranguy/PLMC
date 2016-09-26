module procedures_readers_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use procedures_command_arguments, only: create_filename_from_argument
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

    subroutine create(readers, environment, components)
        type(Readers_Wrapper), intent(out) :: readers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)

        call complete_coordinates_reader_create(readers%complete_coordinates, environment%&
            periodic_box, environment%box_size_checker, components)
    end subroutine create

    subroutine destroy(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call complete_coordinates_reader_destroy(readers%complete_coordinates)
    end subroutine destroy

    subroutine set_from_json(readers, generating_data, prefix)
        type(Readers_Wrapper), intent(in) :: readers
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: coordinates_filename, data_field
        logical :: data_found

        data_field = prefix//"initial coordinates"
        call generating_data%get(data_field, coordinates_filename, data_found)
        call check_data_found(data_field, data_found)
        call readers%complete_coordinates%read(coordinates_filename)
    end subroutine set_from_json

    subroutine set_from_snap(readers, offset_i_snap)
        type(Readers_Wrapper), intent(in) :: readers
        integer, intent(in) :: offset_i_snap

        character(len=:), allocatable :: snap_filename

        call create_filename_from_argument(snap_filename, offset_i_snap)
        call readers%complete_coordinates%read(snap_filename)
    end subroutine set_from_snap

end module procedures_readers_factory
