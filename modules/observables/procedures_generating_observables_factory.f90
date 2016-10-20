module procedures_generating_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
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
        type(Generating_Observables_Wrapper), intent(out) :: observables
        integer, intent(in) :: num_boxes, num_components

        allocate(observables%accessible_domains_size(num_dimensions, num_boxes))
        observables%accessible_domains_size = 0._DP
        call observables_changes_create(observables%teleportations_counters, num_boxes, &
            num_components)
        allocate(observables%teleportations_successes(num_components, num_boxes, num_boxes))
        observables%teleportations_successes = 0._DP
        call observables_changes_create(observables%volumes_change_counter, num_boxes)
        call reals_create(observables%volumes_change_success, num_boxes)

        allocate(observables%nums_particles(num_components, num_boxes))
        observables%nums_particles = 0
        call observables_energies_create(observables%energies, num_boxes, num_components)
        call observables_changes_create(observables%changes, num_boxes, num_components)
    end subroutine create

    pure subroutine destroy(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        call observables_changes_destroy(observables%changes)
        call observables_energies_destroy(observables%energies)
        if (allocated(observables%nums_particles)) deallocate(observables%nums_particles)

        call reals_destroy(observables%volumes_change_success)
        call observables_changes_destroy(observables%volumes_change_counter)
        if (allocated(observables%teleportations_successes)) &
            deallocate(observables%teleportations_successes)
        call observables_changes_destroy(observables%teleportations_counters)
        if (allocated(observables%accessible_domains_size)) &
            deallocate(observables%accessible_domains_size)
    end subroutine destroy

end module procedures_generating_observables_factory
