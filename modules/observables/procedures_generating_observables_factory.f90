module procedures_generating_observables_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Change_Counter, Concrete_Changes_Counter, &
    Concrete_Changes_Success, Concrete_Change_Counters_Line, reset_counters
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_observables_factory_micro, only: create_reals, create_triangle_reals, &
    create_square_reals, destroy_reals, destroy_triangle_reals, destroy_square_reals

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_all
    module procedure :: create_num_particles
    module procedure :: create_changes_counters
    module procedure :: create_triangle_counters
    module procedure :: create_square_counters
    module procedure :: create_changes_successes
    module procedure :: create_reals
    module procedure :: create_triangle_reals
    module procedure :: create_square_reals
end interface create

interface destroy
    module procedure :: destroy_square_reals
    module procedure :: destroy_triangle_reals
    module procedure :: destroy_reals
    module procedure :: destroy_changes_successes
    module procedure :: destroy_square_counters
    module procedure :: destroy_triangle_counters
    module procedure :: destroy_changes_counters
    module procedure :: destroy_num_particles
    module procedure :: destroy_all
end interface destroy

contains

    pure subroutine create_all(observables, num_components)
        type(Generating_Observables_Wrapper), intent(out) :: observables
        integer, intent(in) :: num_components

        call create(observables%num_particles, num_components)
        call create(observables%walls_energies, num_components)
        call create(observables%field_energies, num_components)
        call create(observables%short_energies, num_components)
        call create(observables%dipolar_energies, num_components)
        call create(observables%changes_counters, num_components)
        call create(observables%switches_counters, num_components)
        call create(observables%transmutations_counters, num_components)
        call create(observables%changes_sucesses, num_components)
        call create(observables%switches_successes, num_components)
        call create(observables%transmutations_successes, num_components)
    end subroutine create_all

    pure subroutine destroy_all(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        call destroy(observables%transmutations_successes)
        call destroy(observables%switches_successes)
        call destroy(observables%changes_sucesses)
        call destroy(observables%transmutations_counters)
        call destroy(observables%switches_counters)
        call destroy(observables%changes_counters)
        call destroy(observables%dipolar_energies)
        call destroy(observables%short_energies)
        call destroy(observables%field_energies)
        call destroy(observables%walls_energies)
        call destroy(observables%num_particles)
    end subroutine destroy_all

    pure subroutine create_num_particles(num_particles, num_components)
        integer, allocatable, intent(out) :: num_particles(:)
        integer, intent(in) :: num_components

        allocate(num_particles(num_components))
    end subroutine create_num_particles

    pure subroutine destroy_num_particles(num_particles)
        integer, allocatable, intent(inout) :: num_particles(:)

        if (allocated(num_particles)) deallocate(num_particles)
    end subroutine destroy_num_particles

    pure subroutine create_changes_counters(counters, num_components)
        type(Concrete_Changes_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            call reset_counters(counters(i_counter))
        end do
    end subroutine create_changes_counters

    pure subroutine destroy_changes_counters(counters)
        type(Concrete_Changes_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_changes_counters

    pure subroutine create_triangle_counters(counters, num_components)
        type(Concrete_Change_Counters_Line), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_components

        integer :: i_counter

        allocate(counters(num_components))
        do i_counter = 1, size(counters)
            allocate(counters(i_counter)%line(i_counter))
        end do
        call reset_counters(counters)
    end subroutine create_triangle_counters

    pure subroutine destroy_triangle_counters(counters)
        type(Concrete_Change_Counters_Line), allocatable, intent(inout) :: counters(:)

        integer :: i_counter

        if (allocated(counters)) then
            do i_counter = size(counters), 1, -1
                if (allocated(counters(i_counter)%line)) then
                    deallocate(counters(i_counter)%line)
                end if
            end do
            deallocate(counters)
        end if
    end subroutine destroy_triangle_counters

    pure subroutine create_square_counters(counters, num_components)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:, :)
        integer, intent(in) :: num_components

        allocate(counters(num_components, num_components))
        call reset_counters(counters)
    end subroutine create_square_counters

    pure subroutine destroy_square_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:, :)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_square_counters

    pure subroutine create_changes_successes(successes, num_components)
        type(Concrete_Changes_Success), allocatable, intent(out) :: successes(:)
        integer, intent(in) :: num_components

        allocate(successes(num_components))
    end subroutine create_changes_successes

    pure subroutine destroy_changes_successes(successes)
        type(Concrete_Changes_Success), allocatable, intent(inout) :: successes(:)

        if (allocated(successes)) deallocate(successes)
    end subroutine destroy_changes_successes

end module procedures_generating_observables_factory
