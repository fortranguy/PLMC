module procedures_observables_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_reals_factory, only: reals_create => create, reals_destroy => destroy
use module_changes_success, only: Concrete_Change_Counter, Concrete_Changes_Counter, &
    Concrete_Changes_Success, Concrete_Change_Counter_Line, reset_counters
use types_observables_changes, only: Concrete_Observables_Changes

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_line
    module procedure :: create_element
    module procedure :: create_changes_counters
    module procedure :: create_triangle_counters
    module procedure :: create_square_counters, create_rectangle_counters
    module procedure :: create_teleportations_counters
    module procedure :: create_change_counters
    module procedure :: create_changes_successes
end interface create

interface destroy
    module procedure :: destroy_changes_successes
    module procedure :: destroy_changes_counters
    module procedure :: destroy_teleportations_counters
    module procedure :: destroy_rectangle_counters
    module procedure :: destroy_triangle_counters
    module procedure :: destroy_change_counters
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    pure subroutine create_line(changes, num_boxes, num_components)
        type(Concrete_Observables_Changes), allocatable, intent(out) :: changes(:)
        integer, intent(in) :: num_boxes, num_components

        integer :: i_box

        allocate(changes(num_boxes))
        do i_box = 1, size(changes)
            call create(changes(i_box), num_components)
        end do
    end subroutine create_line

    pure subroutine destroy_line(changes)
        type(Concrete_Observables_Changes), allocatable, intent(inout) :: changes(:)

        integer :: i_box

        if (allocated(changes)) then
            do i_box = size(changes), 1, -1
                call destroy(changes(i_box))
            end do
            deallocate(changes)
        end if
    end subroutine destroy_line

    pure subroutine create_element(changes, num_components)
        type(Concrete_Observables_Changes), intent(out) :: changes
        integer, intent(in) :: num_components

        call create(changes%changes_counters, num_components)
        call create(changes%changes_sucesses, num_components)
        call create(changes%switches_counters, num_components)
        call reals_create(changes%switches_successes, num_components)
        call create(changes%transmutations_counters, num_components)
        allocate(changes%transmutations_successes(num_components, num_components))
        changes%transmutations_successes = 0._DP
    end subroutine create_element

    pure subroutine destroy_element(changes)
        type(Concrete_Observables_Changes), intent(out) :: changes

        if (allocated(changes%transmutations_successes)) &
            deallocate(changes%transmutations_successes)
        call destroy(changes%transmutations_counters)
        call reals_destroy(changes%switches_successes)
        call destroy(changes%switches_counters)
        call destroy(changes%changes_sucesses)
        call destroy(changes%changes_counters)
    end subroutine destroy_element

    pure subroutine create_changes_counters(counters, num_elements)
        type(Concrete_Changes_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_elements

        integer :: i_counter

        allocate(counters(num_elements))
        do i_counter = 1, size(counters)
            call reset_counters(counters(i_counter))
        end do
    end subroutine create_changes_counters

    pure subroutine destroy_changes_counters(counters)
        type(Concrete_Changes_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_changes_counters

    pure subroutine create_triangle_counters(counters, num_elements)
        type(Concrete_Change_Counter_Line), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_elements

        integer :: i_counter

        allocate(counters(num_elements))
        do i_counter = 1, size(counters)
            allocate(counters(i_counter)%line(i_counter))
        end do
        call reset_counters(counters)
    end subroutine create_triangle_counters

    pure subroutine destroy_triangle_counters(counters)
        type(Concrete_Change_Counter_Line), allocatable, intent(inout) :: counters(:)

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

    pure subroutine create_square_counters(counters, num_elements)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:, :)
        integer, intent(in) :: num_elements

        allocate(counters(num_elements, num_elements))
        call reset_counters(counters)
    end subroutine create_square_counters

    pure subroutine create_rectangle_counters(counters, num_elements_1, num_elements_2)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:, :)
        integer, intent(in) :: num_elements_1, num_elements_2

        allocate(counters(num_elements_1, num_elements_2))
        call reset_counters(counters)
    end subroutine create_rectangle_counters

    pure subroutine destroy_rectangle_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:, :)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_rectangle_counters

    pure subroutine create_teleportations_counters(counters, num_boxes, num_components)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:, :, :)
        integer, intent(in) :: num_boxes, num_components

        allocate(counters(num_components, num_boxes, num_boxes))
        call reset_counters(counters)
    end subroutine create_teleportations_counters

    pure subroutine destroy_teleportations_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:, :, :)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_teleportations_counters

    pure subroutine create_change_counters(counters, num_elements)
        type(Concrete_Change_Counter), allocatable, intent(out) :: counters(:)
        integer, intent(in) :: num_elements

        allocate(counters(num_elements))
        call reset_counters(counters)
    end subroutine create_change_counters

    pure subroutine destroy_change_counters(counters)
        type(Concrete_Change_Counter), allocatable, intent(inout) :: counters(:)

        if (allocated(counters)) deallocate(counters)
    end subroutine destroy_change_counters

    pure subroutine create_changes_successes(successes, num_elements)
        type(Concrete_Changes_Success), allocatable, intent(out) :: successes(:)
        integer, intent(in) :: num_elements

        allocate(successes(num_elements))
    end subroutine create_changes_successes

    pure subroutine destroy_changes_successes(successes)
        type(Concrete_Changes_Success), allocatable, intent(inout) :: successes(:)

        if (allocated(successes)) deallocate(successes)
    end subroutine destroy_changes_successes

end module procedures_observables_changes_factory
