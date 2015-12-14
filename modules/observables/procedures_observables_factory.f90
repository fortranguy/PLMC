module procedures_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success, &
    counter_reset => Concrete_Changes_Counter_reset
use types_observables_wrapper, only: Concrete_Components_Energies, Observables_Wrapper

implicit none

private
public :: observables_create, observables_destroy, create_components_energies_nodes, &
    destroy_components_energies_nodes

interface observables_create
    module procedure :: create_all
    module procedure :: create_counters
    module procedure :: create_successes
    module procedure :: create_components_energies
    module procedure :: create_energies
end interface observables_create

interface observables_destroy
    module procedure :: destroy_energies
    module procedure :: destroy_components_energies
    module procedure :: destroy_successes
    module procedure :: destroy_counters
    module procedure :: destroy_all
end interface observables_destroy

contains

    pure subroutine create_all(observables, components)
        type(Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:)

        call observables_create(observables%changes_counters, size(components))
        call observables_create(observables%changes_sucesses, size(components))
        call observables_create(observables%short_energies, size(components))
        call observables_create(observables%walls_energies, size(components))
        call observables_create(observables%long_energies_wo_reci, size(components))
        call observables_create(observables%field_energies, size(components))
    end subroutine create_all

    pure subroutine destroy_all(observables)
        type(Observables_Wrapper), intent(inout) :: observables

        call observables_destroy(observables%field_energies)
        call observables_destroy(observables%long_energies_wo_reci)
        call observables_destroy(observables%walls_energies)
        call observables_destroy(observables%short_energies)
        call observables_destroy(observables%changes_sucesses)
        call observables_destroy(observables%changes_counters)
    end subroutine destroy_all

    pure subroutine create_counters(counters, num_components)
        type(Concrete_Changes_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            call counter_reset(counters(i_counter))
        end do
    end subroutine create_counters

    pure subroutine destroy_counters(counters)
        type(Concrete_Changes_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_counters

    pure subroutine create_successes(successes, num_components)
        type(Concrete_Changes_Success), allocatable, intent(out) :: successes(:)
        integer, intent(in) :: num_components

        allocate(successes(num_components))
    end subroutine create_successes

    pure subroutine destroy_successes(successes)
        type(Concrete_Changes_Success), allocatable, intent(inout) :: successes(:)

        if (allocated(successes)) deallocate(successes)
    end subroutine destroy_successes

    pure subroutine create_components_energies(energies, num_components)
        type(Concrete_Components_Energies), allocatable, intent(out) :: energies(:)
        integer, intent(in) :: num_components

        allocate(energies(num_components))
        call create_components_energies_nodes(energies)
    end subroutine create_components_energies

    pure subroutine destroy_components_energies(energies)
        type(Concrete_Components_Energies), allocatable, intent(inout) :: energies(:)

        if (allocated(energies)) then
            call destroy_components_energies_nodes(energies)
            deallocate(energies)
        end if
    end subroutine destroy_components_energies

    pure subroutine create_components_energies_nodes(energies)
        type(Concrete_Components_Energies), intent(out) :: energies(:)

        integer :: i_component
        do i_component = 1, size(energies)
            allocate(energies(i_component)%with_components(i_component))
            energies(i_component)%with_components = 0._DP
        end do
    end subroutine create_components_energies_nodes

    pure subroutine destroy_components_energies_nodes(energies)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)

        integer :: i_component

        do i_component = size(energies), 1, -1
            deallocate(energies(i_component)%with_components)
        end do
    end subroutine destroy_components_energies_nodes

    pure subroutine create_energies(energies, num_components)
        real(DP), allocatable, intent(out) :: energies(:)
        integer, intent(in) :: num_components

        allocate(energies(num_components))
    end subroutine create_energies

    pure subroutine destroy_energies(energies)
        real(DP), allocatable, intent(inout) :: energies(:)

        if (allocated(energies)) deallocate(energies)
    end subroutine destroy_energies

end module procedures_observables_factory
