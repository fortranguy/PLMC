module procedures_short_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_walls_factory, only: walls_create => create, walls_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only:  Mixture_Wrapper
use procedures_pairs_factory, only: pairs_create => create, pairs_destroy => destroy
use classes_visitable_list, only: Abstract_Visitable_List
use procedures_visitable_list_factory, only: visitable_list_allocate => allocate, &
    visitable_list_deallocate => deallocate
use procedures_cells_factory, only: cells_create => create, cells_destroy => destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper

implicit none

private
public :: short_interactions_create, short_interactions_destroy

contains

    subroutine short_interactions_create(short_interactions, environment, mixture, generating_data,&
        prefix)
        type(Short_Interactions_Wrapper), intent(out) :: short_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        logical :: interact_with_walls, interact
        class(Abstract_Visitable_List), allocatable :: list

        call pairs_create(short_interactions%wall_pairs, interact_with_walls, mixture%&
            wall_min_distances, generating_data, prefix)
        call walls_create(short_interactions%walls_visitor, environment%walls, interact_with_walls)
        call pairs_create(short_interactions%components_pairs, interact, mixture%&
            components_min_distances, generating_data, prefix)
        call pairs_create(short_interactions%components_visitor, environment%periodic_box, interact)
        call cells_create(short_interactions%neighbour_cells, environment%periodic_box, &
            short_interactions%components_pairs, interact)
        call visitable_list_allocate(list, interact, generating_data, prefix)
        call cells_create(short_interactions%visitable_cells, environment%periodic_box, mixture%&
            components, short_interactions%components_pairs, short_interactions%neighbour_cells, &
            list, interact)
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
    end subroutine short_interactions_destroy

end module procedures_short_interactions_factory
