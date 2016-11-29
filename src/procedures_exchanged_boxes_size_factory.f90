module procedures_exchanged_boxes_size_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_exchanged_boxes_size, only: XYZ_Exchanged_Boxes_Size, XY_Exchanged_Boxes_Size, &
    Null_Exchanged_Boxes_Size, Exchanged_Boxes_Size_Line
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

contains

    subroutine create(exchanged_boxes_size, periodic_boxes, tuning_parameters, generating_data, &
        prefix)
        type(Exchanged_Boxes_Size_Line), allocatable, intent(out) :: exchanged_boxes_size(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box, j_box
        type(Concrete_Number_to_String) :: string
        real(DP) :: frequency_ratio, initial_delta
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (size(periodic_boxes) > 1) then
            data_field = prefix//"frequency ratio"
            call generating_data%get(data_field, frequency_ratio, data_found)
            call check_data_found(data_field, data_found)
        end if

        allocate(exchanged_boxes_size(size(periodic_boxes)))
        if (all(periodicity_is_xyz(periodic_boxes))) then
            call allocate_xyz(exchanged_boxes_size)
        else if (all(periodicity_is_xy(periodic_boxes))) then
            call allocate_xy(exchanged_boxes_size)
        else
            call error_exit("procedures_exchanged_boxes_size_factory: create: "//&
                "box periodicity is unknown.")
        end if

        do j_box = 1, size(exchanged_boxes_size)
            do i_box = 1, size(exchanged_boxes_size(j_box)%line)
                if (i_box /= j_box) then
                    data_field = prefix//"Boxes "//string%get(i_box)//string%get(j_box)//&
                    ".initial delta"
                    call generating_data%get(data_field, initial_delta, data_found)
                    call check_data_found(data_field, data_found)
                else
                    initial_delta = 0._DP
                end if
                call exchanged_boxes_size(j_box)%line(i_box)%exchanged%&
                    set(frequency_ratio, initial_delta, tuning_parameters)
            end do
        end do
    end subroutine create

    subroutine allocate_xyz(exchanged_boxes_size)
        type(Exchanged_Boxes_Size_Line), intent(inout) :: exchanged_boxes_size(:)

        integer :: i_box, j_box

        do j_box = 1, size(exchanged_boxes_size)
            allocate(exchanged_boxes_size(j_box)%line(j_box))
            do i_box = 1, size(exchanged_boxes_size(j_box)%line)
                if (i_box /= j_box) then
                    allocate(XYZ_Exchanged_Boxes_Size :: exchanged_boxes_size(j_box)%line(i_box)%&
                        exchanged)
                else
                    allocate(Null_Exchanged_Boxes_Size :: exchanged_boxes_size(j_box)%line(i_box)%&
                        exchanged)
                end if
            end do
        end do
    end subroutine allocate_xyz

    subroutine allocate_xy(exchanged_boxes_size)
        type(Exchanged_Boxes_Size_Line), intent(inout) :: exchanged_boxes_size(:)

        integer :: i_box, j_box

        do j_box = 1, size(exchanged_boxes_size)
            allocate(exchanged_boxes_size(j_box)%line(j_box))
            do i_box = 1, size(exchanged_boxes_size(j_box)%line)
                if (i_box /= j_box) then
                    allocate(XY_Exchanged_Boxes_Size :: exchanged_boxes_size(j_box)%line(i_box)%&
                        exchanged)
                else
                    allocate(Null_Exchanged_Boxes_Size :: exchanged_boxes_size(j_box)%line(i_box)%&
                        exchanged)
                end if
            end do
        end do
    end subroutine allocate_xy

    subroutine destroy(exchanged_boxes_size)
        type(Exchanged_Boxes_Size_Line), allocatable, intent(inout) :: exchanged_boxes_size(:)

        integer :: i_box, j_box

        if (allocated(exchanged_boxes_size)) then
            do j_box = size(exchanged_boxes_size), 1, -1
                if (allocated(exchanged_boxes_size(j_box)%line)) then
                    do i_box = size(exchanged_boxes_size(j_box)%line), 1, -1
                        deallocate(exchanged_boxes_size(j_box)%line(i_box)%exchanged)
                    end do
                    deallocate(exchanged_boxes_size(j_box)%line)
                end if
            end do
            deallocate(exchanged_boxes_size)
        end if
    end subroutine destroy

end module procedures_exchanged_boxes_size_factory
