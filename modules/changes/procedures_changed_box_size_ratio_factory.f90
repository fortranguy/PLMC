module procedures_changed_box_size_ratio_factory

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

    subroutine create(changes_box_size_ratio, periodic_box, can_change, input_data, prefix)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(out) :: changes_box_size_ratio
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: can_change
        type(json_file), optional, intent(inout) :: input_data
        character(len=*), optional, intent(in) :: prefix

        call allocate(changes_box_size_ratio, periodic_box, can_change)
        call construct(changes_box_size_ratio, input_data, prefix)
    end subroutine create

    subroutine allocate(changes_box_size_ratio, periodic_box, can_change)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(out) :: changes_box_size_ratio
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: can_change

        if (can_change) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(XYZ_Changed_Box_Size_Ratio :: changes_box_size_ratio)
            else if (periodicity_is_xy(periodic_box)) then
                allocate(XY_Changed_Box_Size_Ratio :: changes_box_size_ratio)
            else
                call error_exit("procedures_changed_box_size_ratio_factory: allocate: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Changed_Box_Size_Ratio :: changes_box_size_ratio)
        end if
    end subroutine allocate

    subroutine construct(changed_box_size_ratio, input_data, prefix)
        class(Abstract_Changed_Box_Size_Ratio), intent(inout) :: changed_box_size_ratio
        type(json_file), optional, intent(inout) :: input_data
        character(len=*), optional, intent(in) :: prefix

        real(DP) :: initial_delta
        logical :: is_xyz_or_xy
        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (changed_box_size_ratio)
            type is (XYZ_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .true.
            type is (XY_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .true.
            type is (Null_Changed_Box_Size_Ratio)
                is_xyz_or_xy = .false.
            class default
                call error_exit("procedures_changed_box_size_ratio_factory: construct: "//&
                    "changed_box_size: type unknown.")
        end select
        if (is_xyz_or_xy .and. present(input_data) .and. present(input_data)) then
            data_field = prefix//"Small Change.initial delta"
            call input_data%get(data_field, initial_delta, data_found)
            call check_data_found(data_field, data_found)
        else
            initial_delta = 1._DP
        end if
        call changed_box_size_ratio%set(initial_delta)
    end subroutine construct

    subroutine destroy(changes_box_size_ratio)
        class(Abstract_Changed_Box_Size_Ratio), allocatable, intent(inout) :: changes_box_size_ratio

        if (allocated(changes_box_size_ratio)) deallocate(changes_box_size_ratio)
    end subroutine destroy

end module procedures_changed_box_size_ratio_factory
