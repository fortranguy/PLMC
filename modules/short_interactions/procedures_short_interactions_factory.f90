module procedures_short_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only:  Mixture_Wrapper
use classes_half_distribution, only: Abstract_Half_Distribution
use procedures_half_distribution_factory, only: half_distribution_create => create, &
    half_distribution_destroy => destroy
use procedures_hard_contact_factory, only: hard_contact_create => create, hard_contact_destroy => &
    destroy
use procedures_pairs_factory, only: pairs_create => create, pairs_destroy => destroy
use classes_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_list_factory, only: visitable_list_allocate => allocate, &
    visitable_list_deallocate => deallocate
use procedures_cells_factory, only: cells_create => create, cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_property_inquirers, only: measure_pressure

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

        logical :: interact_with_walls, interact, pressure_needed
        class(Abstract_Half_Distribution), allocatable :: half_distribution
        class(Abstract_Visitable_List), allocatable :: list

        if (present(exploring_data) .and. present(volume_change_prefix)) then
            pressure_needed = measure_pressure(exploring_data, volume_change_prefix)
        else
            pressure_needed = .false.
        end if

        call half_distribution_create(half_distribution, pressure_needed, exploring_data, &
            volume_change_prefix)
        call hard_contact_create(short_interactions%hard_contact, environment%periodic_box, &
            half_distribution, pressure_needed)
        call half_distribution_destroy(half_distribution)
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
        call visitable_list_allocate(list, interact, generating_data, prefix)
        call cells_create(short_interactions%visitable_cells, environment%periodic_box, mixture%&
            components, short_interactions%hard_contact, short_interactions%components_pairs, &
            short_interactions%neighbour_cells, list, interact)
        call visitable_list_deallocate(list)
    end subroutine short_interactions_create

    subroutine short_interactions_destroy(short_interactions)
        type(Short_Interactions_Wrapper), intent(inout) :: short_interactions

        call cells_destroy(short_interactions%visitable_cells)
        call cells_destroy(short_interactions%neighbour_cells)
        call pairs_destroy(short_interactions%components_visitor)
        call pairs_destroy(short_interactions%components_pairs)
        call walls_destroy(short_interactions%walls_visitor)
        call pairs_destroy(short_interactions%wall_pairs)
        call hard_contact_destroy(short_interactions%hard_contact)
    end subroutine short_interactions_destroy

end module procedures_short_interactions_factory
