module procedures_changes_success_writer_factory

use classes_number_to_string, only: Concrete_Number_to_String
use types_string_wrapper, only: String_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only:  component_can_translate, component_can_rotate, &
    component_can_exchange
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_changes_success_writer_selector, only: Changes_Success_Writer_Selector
use classes_changes_success_writer, only: Abstract_Changes_Success_Writer, &
    Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use types_changes_success_writer_wrapper, only: Changes_Success_Writer_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_rectangle
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_rectangle
end interface destroy

contains

    subroutine create_rectangle(components_changes, paths, changes_components, components)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(out) :: components_changes(:, :)
        type(String_Wrapper), intent(in) :: paths(:)
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        type(Changes_Success_Writer_Selector) :: selector_i
        type(Concrete_Number_to_String) :: string
        integer :: i_box, i_component

        allocate(components_changes(size(changes_components, 1), size(changes_components, 2)))
        do i_box = 1, size(components_changes, 2)
            do i_component = 1, size(components_changes, 1)
                selector_i%write_translations = &
                    component_can_translate(changes_components(i_component, i_box)%&
                    translated_positions)
                selector_i%write_rotations = &
                    component_can_rotate(changes_components(i_component, i_box)%&
                        rotated_orientations)
                selector_i%write_exchanges = &
                    component_can_exchange(components(i_component, i_box)%chemical_potential)
                call create(components_changes(i_component, i_box)%writer, paths(i_box)%string//&
                    "changes_success_"//string%get(i_component)//".out", selector_i)
            end do
        end do
    end subroutine create_rectangle

     subroutine destroy_rectangle(components_changes)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(inout) :: components_changes(:, :)

        integer :: i_box, i_component

        if (allocated(components_changes)) then
            do i_box = size(components_changes, 2), 1, -1
                do i_component = size(components_changes, 1), 1, -1
                    call destroy(components_changes(i_component, i_box)%writer)
                end do
            end do
            deallocate(components_changes)
        end if
    end subroutine destroy_rectangle

    subroutine create_element(changes, filename, selector)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes
        character(len=*), intent(in) :: filename
        type(Changes_Success_Writer_Selector), intent(in) :: selector

        if (selector%write_translations) then
            allocate(Concrete_Changes_Success_Writer :: changes)
        else
            allocate(Null_Changes_Success_Writer :: changes)
        end if
        call changes%construct(filename, selector)
    end subroutine create_element

    subroutine destroy_element(changes)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes

        if (allocated(changes)) then
            call changes%destroy()
            deallocate(changes)
        end if
    end subroutine destroy_element

end module procedures_changes_success_writer_factory
