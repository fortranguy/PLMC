module procedures_exploring_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
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

    pure subroutine create(observables, num_boxes, num_components)
        type(Exploring_Observables_Wrapper), intent(out) ::observables
        integer, intent(in) :: num_boxes, num_components

        allocate(observables%gemc_inv_pow_activities(num_components, num_boxes))
        observables%gemc_inv_pow_activities = 0._DP
        call observables_energies_create(observables%gemc_energies, num_boxes, num_components)
        allocate(observables%gemc_insertion_successes(num_components, num_boxes))
        observables%gemc_insertion_successes = 0._DP
        call observables_changes_create(observables%gemc_insertion_counters, num_components, &
            num_boxes)
    end subroutine create

    pure subroutine destroy(observables)
        type(Exploring_Observables_Wrapper), intent(inout) ::observables

        call observables_changes_destroy(observables%gemc_insertion_counters)
        if (allocated(observables%gemc_insertion_successes)) &
            deallocate(observables%gemc_insertion_successes)
        call observables_energies_destroy(observables%gemc_energies)
        if (allocated(observables%gemc_inv_pow_activities)) &
            deallocate(observables%gemc_inv_pow_activities)
    end subroutine destroy

end module procedures_exploring_observables_factory
