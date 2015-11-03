module module_change_tuning

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: warning_continue

implicit none

private
public :: set_increase_factor

    type, public :: Concrete_Change_Tuning_Parameters
        real(DP) :: increase_factor
        real(DP) :: increase_factor_max
    end type Concrete_Change_Tuning_Parameters

contains

    subroutine set_increase_factor(context, current_increase_factor, tuning_parameters, &
        max_factor_reached)
        character(len=*), intent(in) :: context
        real(DP), intent(inout) :: current_increase_factor
        type(Concrete_Change_Tuning_Parameters), intent(in) :: tuning_parameters
        logical, intent(out) :: max_factor_reached

        current_increase_factor = tuning_parameters%increase_factor * current_increase_factor
        if (current_increase_factor > tuning_parameters%increase_factor_max) then
            call warning_continue(context//": current_increase_factor is too big.")
            current_increase_factor = tuning_parameters%increase_factor_max
            max_factor_reached = .true.
        else
            max_factor_reached = .false.
        end if
    end subroutine set_increase_factor

end module module_change_tuning
