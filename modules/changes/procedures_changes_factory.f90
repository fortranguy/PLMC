module procedures_changes_factory

use json_module, only: json_file
use class_number_to_string, only: Concrete_Number_to_String
use class_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_component_factory, only: changes_component_create, changes_component_destroy

implicit none

private
public :: changes_create, changes_destroy

interface changes_create
    module procedure :: create_all
end interface changes_create

interface changes_destroy
    module procedure :: destroy_all
end interface changes_destroy

contains

    subroutine create_all(changes, periodic_box, components, input_data, prefix)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        type(Concrete_Number_to_String) :: string
        integer :: i_component

        allocate(changes%components(size(components)))
        do i_component = 1, size(changes%components)
            call changes_component_create(changes%components(i_component), periodic_box, &
                components(i_component), input_data, prefix//"Component "//string%get(i_component)&
                //".")
        end do
    end subroutine create_all

    subroutine destroy_all(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        integer :: i_component

        if (allocated(changes%components)) then
            do i_component = size(changes%components), 1, -1
                call changes_component_destroy(changes%components(i_component))
            end do
            deallocate(changes%components)
        end if
    end subroutine destroy_all

end module procedures_changes_factory
