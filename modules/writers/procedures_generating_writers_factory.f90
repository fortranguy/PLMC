module procedures_generating_writers_factory

use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_complete_coordinates_writer_factory, only: complete_coordinates_writer_create => &
    create, complete_coordinates_writer_destroy => destroy
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_factory, only: set_can_exchange
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use procedures_square_writer_factory, only: square_writer_create => create_transmutations, &
    square_writer_destroy => destroy
use procedures_changes_success_writer_factory, only: changes_success_writer_create => create, &
    changes_success_writer_destroy => destroy
use procedures_energies_writers_factory, only: energies_writers_create => create, &
    energies_writers_destroy => destroy
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper

implicit none

private
public :: create, destroy

contains

    subroutine create(writers, environment, wall_pairs, components, short_pairs, &
        changes_components, generating_data, prefix)
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        logical, dimension(size(components)) :: can_exchange

        call set_can_exchange(can_exchange, components)
        call line_writer_create(writers%nums_particles, can_exchange, "nums_particles.out")
        call complete_coordinates_writer_create(writers%complete_coordinates, environment%&
            periodic_box, components, "coordinates", generating_data, prefix)
        call energies_writers_create(writers%energies, environment%external_field, wall_pairs, &
            components, short_pairs)
        call changes_success_writer_create(writers%components_changes, changes_components, &
            components)
        call triangle_writer_create(writers%switches_successes, components, &
            "switches_successes.out")
        call square_writer_create(writers%transmutations_successes, components, &
            "transmutations_successes.out")
    end subroutine create

    subroutine destroy(writers)
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call square_writer_destroy(writers%transmutations_successes)
        call triangle_writer_destroy(writers%switches_successes)
        call changes_success_writer_destroy(writers%components_changes)
        call energies_writers_destroy(writers%energies)
        call complete_coordinates_writer_destroy(writers%complete_coordinates)
        call line_writer_destroy(writers%nums_particles)
    end subroutine destroy

end module procedures_generating_writers_factory
