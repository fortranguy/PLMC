module procedures_exploring_readers_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use procedures_command_arguments, only: create_filename_from_argument
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_box_size_reader_factory, only: box_size_reader_create => create, &
    box_size_reader_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use procedures_snap_number_factory, only: snap_number_create => create, &
    snap_number_destroy => destroy
use procedures_coordinates_reader_factory, only: coordinates_reader_create => create, &
    coordinates_reader_destroy => destroy
use types_exploring_readers_wrapper, only: Exploring_Readers_Wrapper

implicit none

private
public :: create, destroy, set

contains

    subroutine create(readers, periodic_box, components, num_snaps, num_offset, generating_data, &
        environment_prefix)
        type(Exploring_Readers_Wrapper), intent(out) :: readers
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        integer, intent(in) :: num_snaps, num_offset
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: environment_prefix

        logical :: box_size_changes
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = environment_prefix//"Box.size changes"
        call generating_data%get(data_field, box_size_changes, data_found)
        call check_data_found(data_field, data_found)

        call box_size_reader_create(readers%box, periodic_box, box_size_changes)
        call snap_number_create(readers%snap_number, num_snaps, num_offset)
        call coordinates_reader_create(readers%components, components)
    end subroutine create

    subroutine destroy(readers)
        type(Exploring_Readers_Wrapper), intent(inout) :: readers

        call coordinates_reader_destroy(readers%components)
        call snap_number_destroy(readers%snap_number)
        call box_size_reader_destroy(readers%box)
    end subroutine destroy

    subroutine set(readers, i_snap)
        type(Exploring_Readers_Wrapper), intent(inout) :: readers
        integer, intent(in) :: i_snap

        character(len=:), allocatable :: snap_filename
        integer :: i_component

        !box size ?
        do i_component = 1, size(readers%components)
            call create_filename_from_argument(snap_filename, readers%snap_number%&
                get(i_component, i_snap))
            !call readers%components(i_component)%coordinates%read(snap_filename)
        end do
    end subroutine set

end module procedures_exploring_readers_factory
