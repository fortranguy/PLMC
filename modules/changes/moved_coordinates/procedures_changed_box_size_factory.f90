module procedures_changed_box_size_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_changed_box_size, only: Abstract_Changed_Box_Size, XYZ_Changed_Box_Size, &
    XY_Changed_Box_Size, Null_Changed_Box_Size
use procedures_checks, only: check_data_found
use module_move_tuning, only: Concrete_Move_Tuning_Parameters
use procedures_property_inquirers, only: periodicity_is_xyz, periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(changed_box_size, periodic_box, tuning_parameters, can_change, &
        generating_data, prefix)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_box_size
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(in) :: can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        call allocate(changed_box_size, periodic_box, can_change)
        call construct(changed_box_size, tuning_parameters, generating_data, prefix)
    end subroutine create

    subroutine allocate(changed_box_size, periodic_box, can_change)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_box_size
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: can_change

        if (can_change) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(XYZ_Changed_Box_Size :: changed_box_size)
            else if (periodicity_is_xy(periodic_box)) then
                allocate(XY_Changed_Box_Size :: changed_box_size)
            else
                call error_exit("procedures_changed_box_size_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Changed_Box_Size :: changed_box_size)
        end if
    end subroutine allocate

    subroutine construct(changed_box_size, tuning_parameters, generating_data, prefix)
        class(Abstract_Changed_Box_Size), intent(inout) :: changed_box_size
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP) :: initial_delta
        logical :: is_xyz_or_xy
        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (changed_box_size)
            type is (XYZ_Changed_Box_Size)
                is_xyz_or_xy = .true.
            type is (XY_Changed_Box_Size)
                is_xyz_or_xy = .true.
            type is (Null_Changed_Box_Size)
                is_xyz_or_xy = .false.
            class default
                call error_exit("procedures_changed_box_size_factory: construct: changed_box_size:"&
                    //" type unknown.")
        end select
        if (is_xyz_or_xy) then
            data_field = prefix//"Small Change.initial delta"
            call generating_data%get(data_field, initial_delta, data_found)
            call check_data_found(data_field, data_found)
        end if
        call changed_box_size%construct(initial_delta, tuning_parameters)
    end subroutine construct

    subroutine destroy(changed_box_size)
        class(Abstract_Changed_Box_Size), allocatable, intent(inout) :: changed_box_size

        if (allocated(changed_box_size)) then
            call changed_box_size%destroy()
            deallocate(changed_box_size)
        end if
    end subroutine destroy

end module procedures_changed_box_size_factory
