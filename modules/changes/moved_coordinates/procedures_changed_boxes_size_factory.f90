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
    Null_Changed_Box_Size, Changed_Box_Size_Line
use procedures_environment_inquirers, only: gemc_box_size_can_change
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_triangle
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_triangle
end interface destroy

contains

    !> If the total volume can change, each box can change its size independently,
    !> cf. Panagiotopoulos Mol. Phys. 1988.
    !> Otherwise, we consider only the half hetero couples, e.g. 1<->2.
    subroutine create_triangle(changed_boxes_size, periodic_boxes, tuning_parameters, &
        total_volume_can_change, generating_data, prefix)
        type(Changed_Box_Size_Line), allocatable, intent(out) :: changed_boxes_size(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(in) :: total_volume_can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box, j_box
        class(Abstract_Changed_Box_Size_Ratio), allocatable :: changed_boxes_size_ratio
        logical :: box_size_can_change
        character(len=:), allocatable :: box_name
        type(Concrete_Number_to_String) :: string

        allocate(changed_boxes_size(size(periodic_boxes)))
        do j_box = 1, size(changed_boxes_size)
            allocate(changed_boxes_size(j_box)%line(j_box))
            do i_box = 1, size(changed_boxes_size(j_box)%line)
                if (total_volume_can_change) then
                    box_size_can_change = i_box == j_box
                    box_name = "Box "//string%get(i_box)//"."
                else if (size(periodic_boxes) > 1) then
                    box_size_can_change = i_box < j_box
                    box_name = "Boxes "//string%get(i_box)//string%get(j_box)//"."
                else
                    box_size_can_change = .false.
                    box_name = ""
                end if
                call changed_boxes_size_ratio_create(changed_boxes_size_ratio, periodic_boxes, &
                    box_size_can_change, generating_data, prefix//box_name)
                call create(changed_boxes_size(j_box)%line(i_box)%changed, &
                    changed_boxes_size_ratio, tuning_parameters, generating_data, prefix//box_name)
                call changed_boxes_size_ratio_destroy(changed_boxes_size_ratio)
            end do
        end do
    end subroutine create_triangle

    subroutine create_element(changed_box_size, changed_box_size_ratio, tuning_parameters, &
        generating_data, prefix)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_box_size
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: changed_box_size_ratio
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP) :: frequency_ratio
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (gemc_box_size_can_change(changed_box_size_ratio)) then
            allocate(Concrete_Changed_Box_Size :: changed_box_size)
            data_field = prefix//"frequency ratio"
            call generating_data%get(data_field, frequency_ratio, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Changed_Box_Size :: changed_box_size)
            frequency_ratio = 0._DP
        end if
        call changed_box_size%construct(changed_box_size_ratio, frequency_ratio, tuning_parameters)
    end subroutine create_element

    subroutine destroy_triangle(changed_boxes_size)
        type(Changed_Box_Size_Line), allocatable, intent(inout) :: changed_boxes_size(:)

        integer :: i_box, j_box

        if (allocated(changed_boxes_size)) then
            do j_box = size(changed_boxes_size), 1, -1
                if (allocated(changed_boxes_size(j_box)%line)) then
                    do i_box = size(changed_boxes_size(j_box)%line), 1, -1
                        call destroy(changed_boxes_size(j_box)%line(i_box)%changed)
                    end do
                    deallocate(changed_boxes_size(j_box)%line)
                end if
            end do
            deallocate(changed_boxes_size)
        end if
    end subroutine destroy_triangle

    subroutine destroy_element(changed_box_size)
        class(Abstract_Changed_Box_Size), allocatable, intent(inout) :: changed_box_size

        if (allocated(changed_box_size)) then
            call changed_box_size%destroy()
            deallocate(changed_box_size)
        end if
    end subroutine destroy_element

end module procedures_changed_boxes_size_factory
