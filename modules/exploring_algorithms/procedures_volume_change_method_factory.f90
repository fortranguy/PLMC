module procedures_volume_change_method_factory


use data_input_prefixes, only: volume_change_prefix
use json_module, only: json_file
use procedures_checks, only: check_data_found
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use classes_volume_change_method, only: Abstract_Volume_Change_Method, &
    Concrete_Volume_Change_Method, Null_Volume_Change_Method

implicit none

private
public :: create, destroy

contains

    subroutine create(volume_change_method, physical_model, changed_boxes_size_ratio, &
        measure_pressure_excess, exploring_data)
        class(Abstract_Volume_Change_Method), allocatable, intent(out) :: volume_change_method
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_boxes_size_ratio(:)
        logical, intent(in) :: measure_pressure_excess
        type(json_file), intent(inout) :: exploring_data

        integer :: num_changes
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (measure_pressure_excess) then
            allocate(Concrete_Volume_Change_Method :: volume_change_method)
            data_field = volume_change_prefix//"number"
            call exploring_data%get(data_field, num_changes, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Volume_Change_Method :: volume_change_method)
            num_changes = 0
        end if
        call volume_change_method%construct(physical_model%environment, physical_model%mixture%&
            components, physical_model%short_interactions, physical_model%&
            dipolar_interactions_facades, changed_boxes_size_ratio, num_changes)
    end subroutine create

    subroutine destroy(volume_change_method)
        class(Abstract_Volume_Change_Method), allocatable, intent(inout) :: volume_change_method

        if (allocated(volume_change_method)) then
            call volume_change_method%destroy()
            deallocate(volume_change_method)
        end if
    end subroutine destroy

end module procedures_volume_change_method_factory
