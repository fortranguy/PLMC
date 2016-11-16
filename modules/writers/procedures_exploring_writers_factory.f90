module procedures_exploring_writers_factory

use data_input_prefixes, only: writers_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_pair_potential, only: Pair_Potential_Wrapper, Pair_Potential_Line
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_energies_writers_factory, only: energies_writers_create => create, &
    energies_writers_destroy => destroy
use procedures_directed_graph_writer_factory, only: directed_graph_writer_create => create, &
    directed_graph_writer_destroy => destroy
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure, &
    measure_chemical_potentials, make_dipoles_graph

implicit none

private
public :: create, destroy

contains

    !> @note cf. [[procedures_generating_writers_factory:create]] for a comment about box_stat_i.
    subroutine create(writers, environment, wall_pairs, components, short_pairs, &
        markov_chain_explorer, visit_energies, generating_data)
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(Pair_Potential_Line), intent(in) :: short_pairs(:)
        type(Markov_Chain_Explorer_Wrapper), intent(in) :: markov_chain_explorer
        logical, intent(in) :: visit_energies
        type(json_file), intent(inout) :: generating_data

        type(String_Wrapper), dimension(size(environment%periodic_boxes)) :: boxes_path, graphs_path
        logical :: selectors(size(components, 1), size(components, 2)), measure_pressure
        character(len=:), allocatable :: make_directory_cmd, separator
        logical :: write_dipoles_graph, data_found
        character(len=:), allocatable :: data_field
        type(Concrete_Number_to_String) :: string
        integer :: i_box, box_stat_i

        write_dipoles_graph = all(make_dipoles_graph(markov_chain_explorer%&
            dipolar_neighbourhoods_visitors))
        data_field = writers_prefix//"Shell.make directory command"
        call generating_data%get(data_field, make_directory_cmd, data_found)
        call check_data_found(data_field, data_found)
        data_field = writers_prefix//"Shell.path separator"
        call generating_data%get(data_field, separator, data_found)
        call check_data_found(data_field, data_found)
        box_stat_i = 1
        do i_box = 1, size(boxes_path)
            boxes_path(i_box)%string = "exploring_box_"//string%get(i_box)//separator
            call execute_command_line(make_directory_cmd//" "//boxes_path(i_box)%string, &
                exitstat=box_stat_i)
            if (box_stat_i /= 0) call error_exit("procedures_exploring_writers_factory: create:"//&
                " "//boxes_path(i_box)%string//" directory can't be created.")
            if (.not.write_dipoles_graph) cycle
            graphs_path(i_box)%string = boxes_path(i_box)%string//"graphs"//separator
            call execute_command_line(make_directory_cmd//" "//graphs_path(i_box)%string, &
                exitstat=box_stat_i)
            if (box_stat_i /= 0) call error_exit("procedures_exploring_writers_factory: create:"//&
                " "//graphs_path(i_box)%string//" directory can't be created.")
        end do

        measure_pressure = property_measure_pressure(markov_chain_explorer%volume_change_method)
        call real_writer_create(writers%maximum_boxes_compression_delta, boxes_path, &
            "maximum_box_compression_delta.out", measure_pressure)
        call real_writer_create(writers%beta_pressures_excess, boxes_path, &
            "beta_pressure_excess.out", measure_pressure)
        call energies_writers_create(writers%energies, boxes_path, environment%external_fields, &
            wall_pairs, components, short_pairs, visit_energies)
        selectors = measure_chemical_potentials(markov_chain_explorer%particle_insertion_method)
        call line_writer_create(writers%insertion_successes, boxes_path, "insertion_successes.out",&
            selectors)
        call line_writer_create(writers%inv_pow_activities, boxes_path, "inv_pow_activities.out", &
            selectors)
        call directed_graph_writer_create(writers%dipoles_graph_writer, graphs_path, &
            "dipoles_graph", "dipoles", write_dipoles_graph)
    end subroutine create

    subroutine destroy(writers)
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call directed_graph_writer_destroy(writers%dipoles_graph_writer)
        call energies_writers_destroy(writers%energies)
        call line_writer_destroy(writers%inv_pow_activities)
        call line_writer_destroy(writers%insertion_successes)
        call real_writer_destroy(writers%beta_pressures_excess)
        call real_writer_destroy(writers%maximum_boxes_compression_delta)
    end subroutine destroy

end module procedures_exploring_writers_factory
