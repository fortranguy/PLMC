module procedures_volume_change_method_factory

use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_volume_change_method, only: Abstract_Volume_Change_Method, &
    Concrete_Volume_Change_Method, Null_Volume_Change_Method

implicit none

private
public :: create, destroy

contains

    subroutine create(volume_change_method, physical_model, changed_box_size_ratio, &
        measure_pressure_excess)
        class(Abstract_Volume_Change_Method), allocatable, intent(out) :: volume_change_method
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_box_size_ratio
        logical, intent(in) :: measure_pressure_excess

        if (measure_pressure_excess) then
            allocate(Concrete_Volume_Change_Method :: volume_change_method)
        else
            allocate(Null_Volume_Change_Method :: volume_change_method)
        end if
        call volume_change_method%construct(physical_model%environment, physical_model%mixture%&
            components, physical_model%short_interactions, changed_box_size_ratio)
    end subroutine create

    subroutine destroy(volume_change_method)
        class(Abstract_Volume_Change_Method), allocatable, intent(inout) :: volume_change_method

        if (allocated(volume_change_method)) then
            call volume_change_method%destroy()
            deallocate(volume_change_method)
        end if
    end subroutine destroy

end module procedures_volume_change_method_factory
