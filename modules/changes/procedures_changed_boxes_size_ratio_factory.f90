module procedures_changed_boxes_size_ratio_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_inquirers, only: periodicity_is_xyz, periodicity_is_xy
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio, &
    XYZ_Changed_Box_Size_Ratio, XY_Changed_Box_Size_Ratio, Null_Changed_Box_Size_Ratio

implicit none

private
public :: create, destroy

contains

    subroutine create(changed_boxes_size_ratio, periodic_boxes, volume_can_change, input_data, &
        prefix)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(out) :: &
            changed_boxes_size_ratio(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        logical, intent(in) :: volume_can_change
        type(json_file), optional, intent(inout) :: input_data
        character(len=*), optional, intent(in) :: prefix

        call allocate(changed_boxes_size_ratio, periodic_boxes, volume_can_change)
        call construct(changed_boxes_size_ratio, input_data, prefix)
    end subroutine create

    subroutine allocate(changed_boxes_size_ratio, periodic_boxes, volume_can_change)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(out) :: &
            changed_boxes_size_ratio(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        logical, intent(in) :: volume_can_change

        if (volume_can_change) then
            if (all(periodicity_is_xyz(periodic_boxes))) then
                allocate(XYZ_Changed_Box_Size_Ratio :: &
                    changed_boxes_size_ratio(size(periodic_boxes)))
            else if (all(periodicity_is_xy(periodic_boxes))) then
                allocate(XY_Changed_Box_Size_Ratio :: &
                    changed_boxes_size_ratio(size(periodic_boxes)))
            else
                call error_exit("procedures_changed_boxes_size_ratio_factory: allocate: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Changed_Box_Size_Ratio :: changed_boxes_size_ratio(size(periodic_boxes)))
        end if
    end subroutine allocate

    subroutine construct(changed_boxes_size_ratio, input_data, prefix)
        class(Abstract_Changed_Box_Size_Ratio), intent(inout) :: changed_boxes_size_ratio(:)
        type(json_file), optional, intent(inout) :: input_data
        character(len=*), optional, intent(in) :: prefix

        integer :: i_box
        real(DP) :: initial_delta
        logical :: is_xyz_or_xy
        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (changed_boxes_size_ratio)
            type is (XYZ_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .true.
            type is (XY_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .true.
            type is (Null_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .false.
            class default
                call error_exit("procedures_changed_boxes_size_ratio_factory: construct: "//&
                    "changed_box_size: type unknown.")
        end select
        if (present(input_data) .and. present(prefix) .and. is_xyz_or_xy) then
            data_field = prefix//"Small Change.initial delta"
            call input_data%get(data_field, initial_delta, data_found)
            call check_data_found(data_field, data_found)
        else
            initial_delta = 1._DP
        end if

        do i_box = 1, size(changed_boxes_size_ratio)
            call changed_boxes_size_ratio(i_box)%set(initial_delta)
        end do
    end subroutine construct

    subroutine destroy(changed_boxes_size_ratio)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(inout) :: &
            changed_boxes_size_ratio(:)

        if (allocated(changed_boxes_size_ratio)) deallocate(changed_boxes_size_ratio)
    end subroutine destroy

end module procedures_changed_boxes_size_ratio_factory
