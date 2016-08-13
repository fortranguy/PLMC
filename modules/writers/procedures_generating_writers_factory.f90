module procedures_generating_writers_factory

use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use procedures_complete_coordinates_writer_factory, only: complete_coordinates_writer_create => &
    create, complete_coordinates_writer_destroy => destroy
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_factory, only: set_can_exchange
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use procedures_square_writer_factory, only: square_writer_create => create_transmutations, &
    square_writer_destroy => destroy
use procedures_changes_success_writer_factory, only: changes_success_writer_create => create, &
    changes_success_writer_destroy => destroy
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

        logical, dimension(size(components)) :: can_exchange, are_dipolar

        call line_writer_create(writers%field, environment%external_field, components, &
            "field_energies.out")
        call line_writer_create(writers%walls, wall_pairs, "walls_energies.out")
        call set_can_exchange(can_exchange, components)
        call line_writer_create(writers%num_particles, can_exchange, "num_particles.out")
        call complete_coordinates_writer_create(writers%complete_coordinates, environment%&
            periodic_box, components, "coordinates", generating_data, prefix)
        call changes_success_writer_create(writers%components_changes, changes_components, &
            components)
        call triangle_writer_create(writers%short_energies, short_pairs, "short_energies.out")
        call set_are_dipolar(are_dipolar, components)
        call triangle_writer_create(writers%dipolar_energies, are_dipolar, "dipolar_energies.out")
        call real_writer_create(writers%dipolar_mixture_energy, any(are_dipolar), &
            "dipolar_mixture_energy.out")
        call triangle_writer_create(writers%switches, components, "switches.out")
        call square_writer_create(writers%transmutations, components, "transmutations.out")
    end subroutine create

    subroutine destroy(writers)
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call square_writer_destroy(writers%transmutations)
        call triangle_writer_destroy(writers%switches)
        call real_writer_destroy(writers%dipolar_mixture_energy)
        call triangle_writer_destroy(writers%dipolar_energies)
        call triangle_writer_destroy(writers%short_energies)
        call changes_success_writer_destroy(writers%components_changes)
        call complete_coordinates_writer_destroy(writers%complete_coordinates)
        call line_writer_destroy(writers%num_particles)
        call line_writer_destroy(writers%walls)
        call line_writer_destroy(writers%field)
    end subroutine destroy

end module procedures_generating_writers_factory
