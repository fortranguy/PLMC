module procedures_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use types_reals_line, only: Reals_Line
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success, &
    Concrete_Switch_Counters, changes_counter_reset, switches_counters_reset
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: observables_create, observables_destroy, create_triangle_nodes, destroy_triangle_nodes

interface observables_create
    module procedure :: create_all
    module procedure :: create_changes_counters
    module procedure :: create_switches_counters
    module procedure :: create_changes_successes
    module procedure :: create_triangle
    module procedure :: create_energies
end interface observables_create

interface observables_destroy
    module procedure :: destroy_energies
    module procedure :: destroy_triangle
    module procedure :: destroy_changes_successes
    module procedure :: destroy_switches_counters
    module procedure :: destroy_changes_counters
    module procedure :: destroy_all
end interface observables_destroy

contains

    pure subroutine create_all(observables, components)
        type(Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:)

        call observables_create(observables%changes_counters, size(components))
        call observables_create(observables%switches_counters, size(components))
        call observables_create(observables%changes_sucesses, size(components))
        call observables_create(observables%switches_successes, size(components))
        call observables_create(observables%walls_energies, size(components))
        call observables_create(observables%field_energies, size(components))
        call observables_create(observables%short_energies, size(components))
        call observables_create(observables%dipolar_energies, size(components))
    end subroutine create_all

    pure subroutine destroy_all(observables)
        type(Observables_Wrapper), intent(inout) :: observables

        call observables_destroy(observables%dipolar_energies)
        call observables_destroy(observables%short_energies)
        call observables_destroy(observables%field_energies)
        call observables_destroy(observables%walls_energies)
        call observables_destroy(observables%switches_successes)
        call observables_destroy(observables%changes_sucesses)
        call observables_destroy(observables%switches_counters)
        call observables_destroy(observables%changes_counters)
    end subroutine destroy_all

    pure subroutine create_changes_counters(counters, num_components)
        type(Concrete_Changes_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            call changes_counter_reset(counters(i_counter))
        end do
    end subroutine create_changes_counters

    pure subroutine destroy_changes_counters(counters)
        type(Concrete_Changes_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_changes_counters

    pure subroutine create_switches_counters(counters, num_components)
        type(Concrete_Switch_Counters), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            allocate(counters(i_counter)%line(i_counter))
        end do
        call switches_counters_reset(counters)
    end subroutine create_switches_counters

    pure subroutine destroy_switches_counters(counters)
        type(Concrete_Switch_Counters), allocatable, intent(out) :: counters(:)

        integer :: i_counter

        if (allocated(counters)) then
            do i_counter = size(counters), 1, -1
                if (allocated(counters(i_counter)%line)) then
                    deallocate(counters(i_counter)%line)
                end if
            end do
            deallocate(counters)
        end if
    end subroutine destroy_switches_counters

    pure subroutine create_changes_successes(successes, num_components)
        type(Concrete_Changes_Success), allocatable, intent(out) :: successes(:)
        integer, intent(in) :: num_components

        allocate(successes(num_components))
    end subroutine create_changes_successes

    pure subroutine destroy_changes_successes(successes)
        type(Concrete_Changes_Success), allocatable, intent(inout) :: successes(:)

        if (allocated(successes)) deallocate(successes)
    end subroutine destroy_changes_successes

    pure subroutine create_triangle(triangle, num_components)
        type(Reals_Line), allocatable, intent(out) :: triangle(:)
        integer, intent(in) :: num_components

        allocate(triangle(num_components))
        call create_triangle_nodes(triangle)
    end subroutine create_triangle

    pure subroutine destroy_triangle(triangle)
        type(Reals_Line), allocatable, intent(inout) :: triangle(:)

        if (allocated(triangle)) then
            call destroy_triangle_nodes(triangle)
            deallocate(triangle)
        end if
    end subroutine destroy_triangle

    pure subroutine create_triangle_nodes(triangle)
        type(Reals_Line), intent(out) :: triangle(:)

        integer :: i_component
        do i_component = 1, size(triangle)
            allocate(triangle(i_component)%line(i_component))
            triangle(i_component)%line = 0._DP
        end do
    end subroutine create_triangle_nodes

    pure subroutine destroy_triangle_nodes(triangle)
        type(Reals_Line), intent(inout) :: triangle(:)

        integer :: i_component

        do i_component = size(triangle), 1, -1
            if (allocated(triangle(i_component)%line)) then
                deallocate(triangle(i_component)%line)
            end if
        end do
    end subroutine destroy_triangle_nodes

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
