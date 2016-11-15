module procedures_generating_writers_factory

use data_input_prefixes, only: writers_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_complete_coordinates_writer_factory, only: complete_coordinates_writer_create => &
    create, complete_coordinates_writer_destroy => destroy
use classes_pair_potential, only: Pair_Potential_Wrapper, Pair_Potential_Line
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: set_can_translate, set_can_exchange
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use procedures_rectangle_writer_factory, only: rectangle_writer_create => create, &
    rectangle_writer_destroy => destroy
use procedures_changes_success_writer_factory, only: changes_success_writer_create => create, &
    changes_success_writer_destroy => destroy
use procedures_energies_writers_factory, only: energies_writers_create => create, &
    energies_writers_destroy => destroy
use procedures_writers_inquirers, only: property_write_coordinates => write_coordinates
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper

implicit none

private
public :: create, destroy

contains

    !> @note Valgrind seems to suggest to initialise box_stat_i.
    subroutine create(writers, environment, wall_pairs, components, short_pairs, changes, &
        generating_data)
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Pair_Potential_Line), intent(in) :: short_pairs(:)
        type(Changes_Wrapper), intent(in) :: changes
        type(json_file), intent(inout) :: generating_data

        type(String_Wrapper), dimension(size(environment%periodic_boxes)) :: boxes_path, &
            coordinates_path
        logical :: write_coordinates
        logical :: can_translate(size(changes%components, 1), size(changes%components, 2))
        logical :: can_exchange(size(components, 1), size(components, 2))
        character(len=:), allocatable :: make_directory_cmd, separator
        logical :: data_found
        character(len=:), allocatable :: data_field
        type(Concrete_Number_to_String) :: string
        integer :: i_box, box_stat_i

        write_coordinates = property_write_coordinates(generating_data, writers_prefix)
        data_field = writers_prefix//"Shell.make directory command"
        call generating_data%get(data_field, make_directory_cmd, data_found)
        call check_data_found(data_field, data_found)
        data_field = writers_prefix//"Shell.path separator"
        call generating_data%get(data_field, separator, data_found)
        call check_data_found(data_field, data_found)
        box_stat_i = 1
        do i_box = 1, size(boxes_path)
            boxes_path(i_box)%string = "generating_box_"//string%get(i_box)//separator
            call execute_command_line(make_directory_cmd//" "//boxes_path(i_box)%string, &
                exitstat=box_stat_i)
            if (box_stat_i /= 0) call error_exit("procedures_generating_writers_factory: create:"//&
                " "//boxes_path(i_box)%string//" directory can't be created.")
            if (.not. write_coordinates) cycle
            coordinates_path(i_box)%string = boxes_path(i_box)%string//"coordinates"//separator
            call execute_command_line(make_directory_cmd//" "//coordinates_path(i_box)%string, &
                exitstat=box_stat_i)
            if (box_stat_i /= 0) call error_exit("procedures_generating_writers_factory: create:"//&
                " "//coordinates_path(i_box)%string//" directory can't be created.")
        end do

        call real_writer_create(writers%accessible_domains_size, boxes_path, &
            "accessible_domain_size.out", changes%changed_boxes_size)
        call set_can_translate(can_translate, changes%components)
        call line_writer_create(writers%volumes_change_success, "volumes_change_success.out", &
            changes%changed_boxes_size)
        call triangle_writer_create(writers%volumes_exchange_success, &
            "volumes_exchange_success.out", changes%exchanged_boxes_size)
        call line_writer_create(writers%teleportations_successes, make_directory_cmd, separator, &
            "teleportations", can_translate)
        call rectangle_writer_create(writers%swaps_successes, make_directory_cmd, separator, &
            "swaps", can_translate)

        call set_can_exchange(can_exchange, components)
        call line_writer_create(writers%nums_particles, boxes_path, "nums_particles.out", &
            can_exchange .or. size(can_exchange, 2) > 1)
        call complete_coordinates_writer_create(writers%complete_coordinates, coordinates_path, &
            "coordinates", environment%periodic_boxes, components, generating_data, writers_prefix)

        call energies_writers_create(writers%energies, boxes_path, environment%external_fields, &
            wall_pairs, components, short_pairs, visit_energies=.true.)
        call changes_success_writer_create(writers%components_changes, boxes_path, &
            changes%components, components)
        call triangle_writer_create(writers%switches_successes, boxes_path, &
            "switches_successes.out", components)
        call rectangle_writer_create(writers%transmutations_successes, boxes_path, &
            "transmutations_successes.out", components)
    end subroutine create

    subroutine destroy(writers)
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call rectangle_writer_destroy(writers%transmutations_successes)
        call triangle_writer_destroy(writers%switches_successes)
        call changes_success_writer_destroy(writers%components_changes)

        call energies_writers_destroy(writers%energies)

        call complete_coordinates_writer_destroy(writers%complete_coordinates)
        call line_writer_destroy(writers%nums_particles)

        call rectangle_writer_destroy(writers%swaps_successes)
        call line_writer_destroy(writers%teleportations_successes)
        call triangle_writer_destroy(writers%volumes_exchange_success)
        call line_writer_destroy(writers%volumes_change_success)
        call real_writer_destroy(writers%accessible_domains_size)
    end subroutine destroy

end module procedures_generating_writers_factory
