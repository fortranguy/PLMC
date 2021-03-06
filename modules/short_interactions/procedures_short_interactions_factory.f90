module procedures_short_interactions_factory

use data_input_prefixes, only: short_interactions_prefix, volume_change_prefix, dipoles_graph_prefix
use json_module, only: json_file
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_beta_pressures_excess_factory, only: beta_pressures_excess_create => create, &
    beta_pressures_excess_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: total_volume_can_change
use types_mixture_wrapper, only:  Mixture_Wrapper
use classes_dirac_distribution_plus, only: Abstract_Dirac_Distribution_Plus
use procedures_dirac_distribution_plus_factory, only: dirac_distribution_plus_create => create, &
    dirac_distribution_plus_destroy => destroy
use procedures_hard_contact_factory, only: hard_contact_create => create, hard_contact_destroy => &
    destroy
use procedures_pairs_factory, only: pairs_create => create, pairs_destroy => destroy
use procedures_dipolar_neighbourhoods_factory, only: dipolar_neighbourhoods_create => create, &
    dipolar_neighbourhoods_destroy => destroy
use classes_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_list_factory, only: visitable_list_allocate => allocate, &
    visitable_list_deallocate => deallocate
use procedures_cells_factory, only: cells_create => create, cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure, &
    property_make_dipoles_graph => make_dipoles_graph

implicit none

private
public :: create, destroy

contains

    subroutine create(short_interactions, environment, mixture, generating_data, exploring_data)
        type(Short_Interactions_Wrapper), intent(out) :: short_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: generating_data
        type(json_file), optional, intent(inout) :: exploring_data

        logical :: boxes_size_can_change, measure_pressure, make_dipoles_graph
        logical :: interact_with_walls, interact
        class(Abstract_Dirac_Distribution_Plus), allocatable :: dirac_plus
        class(Abstract_Visitable_List), allocatable :: list_mold

        if (present(exploring_data)) then
            measure_pressure = property_measure_pressure(exploring_data, volume_change_prefix)
            make_dipoles_graph = property_make_dipoles_graph(exploring_data, dipoles_graph_prefix)
        else
            measure_pressure = .false.
            make_dipoles_graph = .false.
        end if

        call beta_pressures_excess_create(short_interactions%beta_pressures_excess, environment%&
            periodic_boxes, environment%accessible_domains, measure_pressure)
        call dirac_distribution_plus_create(dirac_plus, measure_pressure, exploring_data, &
            volume_change_prefix)
        call hard_contact_create(short_interactions%hard_contact, environment%periodic_boxes, &
            dirac_plus, measure_pressure)
        call dirac_distribution_plus_destroy(dirac_plus)
        call pairs_create(short_interactions%wall_pairs, interact_with_walls, mixture%&
            wall_min_distances, generating_data, short_interactions_prefix)
        call walls_create(short_interactions%walls_visitors, environment%visitable_walls, &
            interact_with_walls)
        call pairs_create(short_interactions%components_pairs, interact, mixture%&
            components_min_distances, generating_data, short_interactions_prefix)
        call pairs_create(short_interactions%components_visitors, environment%periodic_boxes, &
            interact)
        call dipolar_neighbourhoods_create(short_interactions%dipolar_neighbourhoods, mixture%&
            components, make_dipoles_graph, exploring_data, dipoles_graph_prefix)
        call visitable_list_allocate(list_mold, interact, generating_data, &
            short_interactions_prefix)
        call cells_create(short_interactions%cells, environment%periodic_boxes, environment%&
            accessible_domains, mixture%components, short_interactions%components_pairs, &
            short_interactions%hard_contact, short_interactions%dipolar_neighbourhoods, &
            list_mold, interact)
        boxes_size_can_change = total_volume_can_change(environment%beta_pressure) .or. &
            size(environment%periodic_boxes) > 1
        call cells_create(short_interactions%visitable_cells_memento, list_mold, &
            (boxes_size_can_change .or. measure_pressure) .and. interact)
        call visitable_list_deallocate(list_mold)
    end subroutine create

    subroutine destroy(short_interactions)
        type(Short_Interactions_Wrapper), intent(inout) :: short_interactions

        call cells_destroy(short_interactions%visitable_cells_memento)
        call cells_destroy(short_interactions%cells)
        call dipolar_neighbourhoods_destroy(short_interactions%dipolar_neighbourhoods)
        call pairs_destroy(short_interactions%components_visitors)
        call pairs_destroy(short_interactions%components_pairs)
        call walls_destroy(short_interactions%walls_visitors)
        call pairs_destroy(short_interactions%wall_pairs)
        call hard_contact_destroy(short_interactions%hard_contact)
        call beta_pressures_excess_destroy(short_interactions%beta_pressures_excess)
    end subroutine destroy

end module procedures_short_interactions_factory
