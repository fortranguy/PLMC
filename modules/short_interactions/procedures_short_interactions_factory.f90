module procedures_short_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use procedures_beta_pressure_excess_factory, only: beta_pressure_excess_create => create, &
    beta_pressure_excess_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: box_size_can_change
use types_mixture_wrapper, only:  Mixture_Wrapper
use classes_dirac_distribution_plus, only: Abstract_Dirac_Distribution_Plus
use procedures_dirac_distribution_plus_factory, only: dirac_distribution_plus_create => create, &
    dirac_distribution_plus_destroy => destroy
use procedures_hard_contact_factory, only: hard_contact_create => create, hard_contact_destroy => &
    destroy
use procedures_pairs_factory, only: pairs_create => create, pairs_destroy => destroy
use classes_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_list_factory, only: visitable_list_allocate => allocate, &
    visitable_list_deallocate => deallocate
use procedures_cells_factory, only: cells_create => create, cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure, &
    property_measure_maximum_compression => measure_maximum_compression

implicit none

private
public :: short_interactions_create, short_interactions_destroy

contains

    subroutine short_interactions_create(short_interactions, environment, mixture, generating_data,&
        prefix, exploring_data, volume_change_prefix)
        type(Short_Interactions_Wrapper), intent(out) :: short_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix
        type(json_file), optional, intent(inout) :: exploring_data
        character(len=*), optional, intent(in) :: volume_change_prefix

        logical :: interact_with_walls, interact, measure_maximum_compression, measure_pressure
        class(Abstract_Dirac_Distribution_Plus), allocatable :: dirac_plus
        class(Abstract_Visitable_List), allocatable :: list_mold

        if (present(exploring_data) .and. present(volume_change_prefix)) then
            measure_pressure = property_measure_pressure(exploring_data, volume_change_prefix)
            measure_maximum_compression = property_measure_maximum_compression(exploring_data, &
                volume_change_prefix)
        else
            measure_pressure = .false.
            measure_maximum_compression = .false.
        end if

        call beta_pressure_excess_create(short_interactions%beta_pressure_excess, environment%&
            periodic_box, environment%accessible_domain, measure_pressure)
        call dirac_distribution_plus_create(dirac_plus, measure_pressure, exploring_data, &
            volume_change_prefix)
        call hard_contact_create(short_interactions%hard_contact, environment%periodic_box, &
            dirac_plus, measure_pressure .or. measure_maximum_compression)
        call dirac_distribution_plus_destroy(dirac_plus)
        call pairs_create(short_interactions%wall_pairs, interact_with_walls, mixture%&
            wall_min_distances, generating_data, prefix)
        call walls_create(short_interactions%walls_visitor, environment%visitable_walls, &
            interact_with_walls)
        call pairs_create(short_interactions%components_pairs, interact, mixture%&
            components_min_distances, generating_data, prefix)
        call pairs_create(short_interactions%components_visitor, environment%periodic_box, interact)
        call cells_create(short_interactions%neighbour_cells, environment%periodic_box, &
            environment%accessible_domain, short_interactions%hard_contact, short_interactions%&
            components_pairs, interact)
        call visitable_list_allocate(list_mold, interact, generating_data, prefix)
        call cells_create(short_interactions%visitable_cells, environment%periodic_box, mixture%&
            components, short_interactions%hard_contact, short_interactions%components_pairs, &
            short_interactions%neighbour_cells, list_mold, interact)
        call cells_create(short_interactions%visitable_cells_memento, list_mold, &
            (box_size_can_change(environment%beta_pressure) .or. measure_pressure) .and. interact)
        call visitable_list_deallocate(list_mold)
    end subroutine short_interactions_create

    subroutine short_interactions_destroy(short_interactions)
        type(Short_Interactions_Wrapper), intent(inout) :: short_interactions

        call cells_destroy(short_interactions%visitable_cells_memento)
        call cells_destroy(short_interactions%visitable_cells)
        call cells_destroy(short_interactions%neighbour_cells)
        call pairs_destroy(short_interactions%components_visitor)
        call pairs_destroy(short_interactions%components_pairs)
        call walls_destroy(short_interactions%walls_visitor)
        call pairs_destroy(short_interactions%wall_pairs)
        call hard_contact_destroy(short_interactions%hard_contact)
        call beta_pressure_excess_destroy(short_interactions%beta_pressure_excess)
    end subroutine short_interactions_destroy

end module procedures_short_interactions_factory
