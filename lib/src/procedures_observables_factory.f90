module procedures_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_mixture_wrapper, only: Mixture_Wrapper
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success, &
    counter_reset => Concrete_Changes_Counter_reset
use types_observables_wrapper, only: Concrete_Energies, Observables_Wrapper

implicit none

private

interface observables_create
    module procedure :: create_all
    module procedure :: create_counters
    module procedure :: create_successes
    module procedure :: create_inter_energies
    module procedure :: create_energies
end interface observables_create

interface observables_destroy
    module procedure :: destroy_energies
    module procedure :: destroy_inter_energies
    module procedure :: destroy_successes
    module procedure :: destroy_counters
    module procedure :: destroy_all
end interface observables_destroy

contains

    subroutine create_all(observables, mixture)
        type(Observables_Wrapper), intent(out) :: observables
        type(Mixture_Wrapper), intent(in) :: mixture

        call observables_create(observables%changes_counters, size(mixture%components))
        call observables_create(observables%changes_sucesses, size(mixture%components))
        call observables_create(observables%short_energies, size(mixture%components))
        call observables_create(observables%walls_energies, size(mixture%components))
        call observables_create(observables%long_energies, size(mixture%components))
        call observables_create(observables%field_energies, size(mixture%components))
    end subroutine create_all

    subroutine destroy_all(observables)
        type(Observables_Wrapper), intent(inout) :: observables

        call observables_destroy(observables%field_energies)
        call observables_destroy(observables%long_energies)
        call observables_destroy(observables%walls_energies)
        call observables_destroy(observables%short_energies)
        call observables_destroy(observables%changes_sucesses)
        call observables_destroy(observables%changes_counters)
    end subroutine destroy_all

    subroutine create_counters(counters, num_components)
        type(Concrete_Changes_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            call counter_reset(counters(i_counter))
        end do
    end subroutine create_counters

    subroutine destroy_counters(counters)
        type(Concrete_Changes_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_counters

    subroutine create_successes(successes, num_components)
        type(Concrete_Changes_Success), allocatable, intent(out) :: successes(:)
        integer, intent(in) :: num_components

        allocate(successes(num_components))
    end subroutine create_successes

    subroutine destroy_successes(successes)
        type(Concrete_Changes_Success), allocatable, intent(inout) :: successes(:)

        if (allocated(successes)) deallocate(successes)
    end subroutine destroy_successes

    subroutine create_inter_energies(inter_energies, num_components)
        type(Concrete_Energies), allocatable, intent(out) :: inter_energies(:)
        integer, intent(in) :: num_components

        integer :: i_component

        allocate(inter_energies(num_components))
        do i_component = 1, size(inter_energies)
            allocate(inter_energies(i_component)%with_components(i_component))
            inter_energies(i_component)%with_components = 0._DP
        end do
    end subroutine create_inter_energies

    subroutine destroy_inter_energies(inter_energies)
        type(Concrete_Energies), allocatable, intent(inout) :: inter_energies(:)

        integer :: i_component

        do i_component = size(inter_energies), 1, -1
            deallocate(inter_energies(i_component)%with_components)
        end do
        deallocate(inter_energies)
    end subroutine destroy_inter_energies

    subroutine create_energies(energies, num_components)
        real(DP), allocatable, intent(out) :: energies(:)
        integer, intent(in) :: num_components

        allocate(energies(num_components))
    end subroutine create_energies

    subroutine destroy_energies(energies)
        real(DP), allocatable, intent(inout) :: energies(:)

        if (allocated(energies)) deallocate(energies)
    end subroutine destroy_energies

end module procedures_observables_factory
