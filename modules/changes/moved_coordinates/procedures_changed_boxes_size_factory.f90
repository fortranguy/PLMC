module procedures_changed_boxes_size_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use classes_changed_box_size, only: Abstract_Changed_Box_Size, Concrete_Changed_Box_Size, &
    Null_Changed_Box_Size
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

contains

    subroutine create(changed_boxes_size, changed_boxes_size_ratio, tuning_parameters, &
        volume_can_change, generating_data, prefix)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_boxes_size(:)
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_boxes_size_ratio(:)
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(in) :: volume_can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        integer :: num_changes
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (volume_can_change) then
            allocate(Concrete_Changed_Box_Size :: &
                changed_boxes_size(size(changed_boxes_size_ratio)))
            data_field = prefix//"number"
            call generating_data%get(data_field, num_changes, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Changed_Box_Size :: changed_boxes_size(size(changed_boxes_size_ratio)))
            num_changes = 0
        end if

        do i_box = 1, size(changed_boxes_size)
            call changed_boxes_size(i_box)%construct(changed_boxes_size_ratio(i_box), num_changes, &
                tuning_parameters)
        end do
    end subroutine create

    subroutine destroy(changed_boxes_size)
        class(Abstract_Changed_Box_Size), allocatable, intent(inout) :: changed_boxes_size(:)

        integer :: i_box

        if (allocated(changed_boxes_size)) then
            do i_box = size(changed_boxes_size), 1, -1
                call changed_boxes_size(i_box)%destroy()
            end do
            deallocate(changed_boxes_size)
        end if
    end subroutine destroy

end module procedures_changed_boxes_size_factory
