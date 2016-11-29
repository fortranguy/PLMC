module procedures_exploring_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_logical_factory, only: logical_destroy => destroy
use module_changes_success, only: Concrete_Change_Counter, reset_counters
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use procedures_reals_factory, only:reals_create => create, reals_destroy => destroy
use procedures_observables_changes_factory, only: observables_changes_create => create, &
    observables_changes_destroy => destroy
use procedures_observables_energies_factory, only: observables_energies_create => create, &
    observables_energies_destroy => destroy

implicit none

private
public :: create, destroy

contains

    !> @note adjacency_matrices%rectangle will be allocated in
    !> [[classes_dipolar_neighbourhoods_visitor:Abstract_try]]
    pure subroutine create(observables, num_boxes, num_components)
        type(Exploring_Observables_Wrapper), intent(out) ::observables
        integer, intent(in) :: num_boxes, num_components

        allocate(observables%maximum_boxes_compression_delta(num_boxes))
        observables%maximum_boxes_compression_delta = 0._DP
        allocate(observables%beta_pressures_excess(num_boxes))
        observables%beta_pressures_excess = 0._DP
        allocate(observables%inv_pow_activities(num_components, num_boxes))
        observables%inv_pow_activities = 0._DP
        call observables_energies_create(observables%energies, num_boxes, num_components)
        call observables_changes_create(observables%insertion_counters, num_components, &
            num_boxes)
        allocate(observables%insertion_successes(num_components, num_boxes))
        observables%insertion_successes = 0._DP
        allocate(observables%adjacency_matrices(num_components, num_components, num_boxes))
    end subroutine create

    pure subroutine destroy(observables)
        type(Exploring_Observables_Wrapper), intent(inout) ::observables

        call logical_destroy(observables%adjacency_matrices)
        if (allocated(observables%insertion_successes)) &
            deallocate(observables%insertion_successes)
        call observables_changes_destroy(observables%insertion_counters)
        call observables_energies_destroy(observables%energies)
        if (allocated(observables%inv_pow_activities)) &
            deallocate(observables%inv_pow_activities)
        if (allocated(observables%beta_pressures_excess)) &
            deallocate(observables%beta_pressures_excess)
        if (allocated(observables%maximum_boxes_compression_delta)) &
            deallocate(observables%maximum_boxes_compression_delta)
    end subroutine destroy

end module procedures_exploring_observables_factory
