module procedures_changes_success_writer_factory

use classes_number_to_string, only: Concrete_Number_to_String
use types_component_wrapper, only: Component_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use classes_changes_success_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use types_changes_success_writer_wrapper, only: Changes_Success_Writer_Wrapper
use procedures_property_inquirers, only:  component_can_translate, component_can_rotate, &
    component_can_exchange

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_line(components_changes, changes_components, components)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(out) :: components_changes(:)
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)
        type(Component_Wrapper), intent(in) :: components(:)

        type(Concrete_Changes_Selector) :: selector_i
        type(Concrete_Number_to_String) :: string
        integer :: i_component

        allocate(components_changes(size(changes_components)))
        do i_component = 1, size(components_changes)
            selector_i%write_translations = &
                component_can_translate(changes_components(i_component)%translated_positions)
            selector_i%write_rotations = component_can_rotate(changes_components(i_component)%&
                rotated_orientations)
            selector_i%write_exchanges = component_can_exchange(components(i_component)%&
                chemical_potential)
            call create(components_changes(i_component)%writer, selector_i, "changes_success_"//&
                string%get(i_component)//".out")
        end do
    end subroutine create_line

     subroutine destroy_line(components_changes)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(inout) :: components_changes(:)

        integer :: i_component

        if (allocated(components_changes)) then
            do i_component = size(components_changes), 1, -1
                call destroy(components_changes(i_component)%writer)
            end do
            deallocate(components_changes)
        end if
    end subroutine destroy_line

    subroutine create_element(changes, selector, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes
        type(Concrete_Changes_Selector), intent(in) :: selector
        character(len=*), intent(in) :: filename

        if (selector%write_translations) then
            allocate(Concrete_Changes_Success_Writer :: changes)
        else
            allocate(Null_Changes_Success_Writer :: changes)
        end if
        call changes%construct(selector, filename)
    end subroutine create_element

    subroutine destroy_element(changes)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes

        if (allocated(changes)) then
            call changes%destroy()
            deallocate(changes)
        end if
    end subroutine destroy_element

end module procedures_changes_success_writer_factory
