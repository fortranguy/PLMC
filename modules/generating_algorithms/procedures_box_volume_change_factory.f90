module procedures_box_volume_change_factory

use classes_dipolar_visitor, only: Abstract_Dipolar_Visitor
use procedures_dipolar_visitor_factory, only: dipolar_visitor_create => create, &
    dipolar_visitor_destroy => destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use classes_box_volume_change, only: Abstract_Box_Volume_Change, Concrete_Box_Volume_Change, &
    Null_Box_Volume_Change
use procedures_environment_inquirers, only: property_box_size_can_change => box_size_can_change

implicit none

private
public :: create, destroy

contains

    subroutine create(box_volume_change, physical_model, changed_box_size)
        class(Abstract_Box_Volume_Change), allocatable, intent(out) :: box_volume_change
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Changed_Box_Size), intent(in) :: changed_box_size

        class(Abstract_Dipolar_Visitor), allocatable :: dipolar_visitor
        logical :: box_size_can_change

        box_size_can_change = property_box_size_can_change(changed_box_size)
        if (box_size_can_change) then
            allocate(Concrete_Box_Volume_Change :: box_volume_change)
        else
            allocate(Null_Box_Volume_Change :: box_volume_change)
        end if
        call dipolar_visitor_create(dipolar_visitor, physical_model%environment%periodic_box, &
            box_size_can_change, physical_model%mixture%components, physical_model%&
            dipolar_interactions_dynamic)
        call box_volume_change%construct(physical_model%environment, physical_model%mixture, &
            physical_model%short_interactions, physical_model%dipolar_interactions_static, &
            dipolar_visitor, changed_box_size)
        call dipolar_visitor_destroy(dipolar_visitor)
    end subroutine create

    subroutine destroy(box_volume_change)
        class(Abstract_Box_Volume_Change), allocatable, intent(inout) :: box_volume_change

        if (allocated(box_volume_change)) then
            call box_volume_change%destroy()
            deallocate(box_volume_change)
        end if
    end subroutine destroy

end module procedures_box_volume_change_factory
