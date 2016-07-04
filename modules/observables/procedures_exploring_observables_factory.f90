module procedures_exploring_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Change_Counter, reset_counters
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use procedures_observables_factory_micro, only: create_reals, create_triangle_reals, destroy_reals,&
    destroy_triangle_reals

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_all
    module procedure :: create_triangle_reals
    module procedure :: create_reals
    module procedure :: create_change_counters
end interface create

interface destroy
    module procedure :: destroy_change_counters
    module procedure :: destroy_reals
    module procedure :: destroy_triangle_reals
    module procedure :: destroy_all
end interface destroy

contains

    pure subroutine create_all(observables, num_components)
        type(Exploring_Observables_Wrapper), intent(out) ::observables
        integer, intent(in) :: num_components

        call create(observables%inv_pow_activities, num_components)
        call create(observables%walls_energies, num_components)
        call create(observables%field_energies, num_components)
        call create(observables%short_energies, num_components)
        call create(observables%dipolar_energies, num_components)
        call create(observables%insertion_successes, num_components)
        call create(observables%insertion_counters, num_components)
    end subroutine create_all

    pure subroutine destroy_all(observables)
        type(Exploring_Observables_Wrapper), intent(inout) ::observables

        call destroy(observables%insertion_counters)
        call destroy(observables%insertion_successes)
        call destroy(observables%dipolar_energies)
        call destroy(observables%short_energies)
        call destroy(observables%field_energies)
        call destroy(observables%walls_energies)
        call destroy(observables%inv_pow_activities)
    end subroutine destroy_all

    pure subroutine create_change_counters(counters, num_components)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        allocate(counters(num_components))
        call reset_counters(counters)
    end subroutine create_change_counters

    pure subroutine destroy_change_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_change_counters

end module procedures_exploring_observables_factory
