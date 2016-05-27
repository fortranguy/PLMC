module procedures_exploring_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use module_changes_success, only: Concrete_Change_Counter, reset_counters
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_all
    module procedure :: create_reals
    module procedure :: create_change_counters
end interface create

interface destroy
    module procedure :: destroy_change_counters
    module procedure :: destroy_reals
    module procedure :: destroy_all
end interface destroy

contains

    subroutine create_all(observables, components)
        type(Exploring_Observables_Wrapper), intent(out) ::observables
        type(Component_Wrapper), intent(in) :: components(:)

        call create(observables%inv_pow_activities, size(components))
        call create(observables%widom_successes, size(components))
        call create(observables%widom_counters, size(components))
    end subroutine create_all

    subroutine destroy_all(observables)
        type(Exploring_Observables_Wrapper), intent(inout) ::observables

        call destroy(observables%widom_counters)
        call destroy(observables%widom_successes)
        call destroy(observables%inv_pow_activities)
    end subroutine destroy_all

    subroutine create_reals(reals, num_components)
        real(DP), allocatable, intent(out) :: reals(:)
        integer, intent(in) :: num_components

        allocate(reals(num_components))
    end subroutine create_reals

    subroutine destroy_reals(reals)
        real(DP), allocatable, intent(inout) :: reals(:)

        if (allocated(reals)) deallocate(reals)
    end subroutine destroy_reals

    subroutine create_change_counters(counters, num_components)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        allocate(counters(num_components))
        call reset_counters(counters)
    end subroutine create_change_counters

    subroutine destroy_change_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_change_counters

end module procedures_exploring_observables_factory
