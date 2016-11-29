module procedures_changed_boxes_size_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_changed_boxes_size_ratio_factory, only: changed_boxes_size_ratio_create => create, &
    changed_boxes_size_ratio_destroy => destroy
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use classes_changed_box_size, only: Abstract_Changed_Box_Size, Concrete_Changed_Box_Size, &
    Null_Changed_Box_Size
use procedures_environment_inquirers, only: box_size_can_change
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

contains

    !> If total volume can change, each box can change its size independently in GEMC
    !> cf. Panagiotopoulos Mol. Phys. 1988.
    subroutine create(changed_boxes_size, periodic_boxes, tuning_parameters, &
        total_volume_can_change, generating_data, prefix)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_boxes_size(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(in) :: total_volume_can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Changed_Box_Size_Ratio), allocatable :: changed_box_size_ratio
        integer :: i_box
        type(Concrete_Number_to_String) :: string
        real(DP) :: frequency_ratio
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (total_volume_can_change) then
            allocate(Concrete_Changed_Box_Size :: changed_boxes_size(size(periodic_boxes)))
            data_field = prefix//"frequency ratio"
            call generating_data%get(data_field, frequency_ratio, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Changed_Box_Size :: changed_boxes_size(size(periodic_boxes)))
            frequency_ratio = 0._DP
        end if

        do i_box = 1, size(changed_boxes_size)
            call changed_boxes_size_ratio_create(changed_box_size_ratio, periodic_boxes(i_box), &
                total_volume_can_change, generating_data, prefix//"Box "//string%get(i_box)//".")
            call changed_boxes_size(i_box)%construct(frequency_ratio, changed_box_size_ratio, &
                tuning_parameters)
            call changed_boxes_size_ratio_destroy(changed_box_size_ratio)
        end do
    end subroutine create

    subroutine destroy(changed_boxes_size)
        class(Abstract_Changed_Box_Size), allocatable, intent(inout) :: changed_boxes_size(:)

        integer :: i_box

        if(allocated(changed_boxes_size)) then
            do i_box = size(changed_boxes_size), 1, -1
                call changed_boxes_size(i_box)%destroy()
            end do
            deallocate(changed_boxes_size)
        end if
    end subroutine destroy

end module procedures_changed_boxes_size_factory
