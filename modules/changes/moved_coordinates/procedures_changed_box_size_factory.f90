module procedures_changed_box_size_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_changed_box_size_ratio, only: Abstract_Changed_Box_Size_Ratio
use classes_changed_box_size, only: Abstract_Changed_Box_Size, Concrete_Changed_Box_Size, &
    Null_Changed_Box_Size
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

contains

    subroutine create(changed_box_size, periodic_box, ratio, tuning_parameters, can_change, &
        generating_data, prefix)
        class(Abstract_Changed_Box_Size), allocatable, intent(out) :: changed_box_size
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: ratio
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(in) :: can_change
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        if (can_change) then
            allocate(Concrete_Changed_Box_Size :: changed_box_size)
        else
            allocate(Null_Changed_Box_Size :: changed_box_size)
        end if
        call changed_box_size%construct(ratio, tuning_parameters)
    end subroutine create

    subroutine destroy(changed_box_size)
        class(Abstract_Changed_Box_Size), allocatable, intent(inout) :: changed_box_size

        if (allocated(changed_box_size)) then
            call changed_box_size%destroy()
            deallocate(changed_box_size)
        end if
    end subroutine destroy

end module procedures_changed_box_size_factory
