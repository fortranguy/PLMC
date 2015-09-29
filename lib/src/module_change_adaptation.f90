module module_change_adaptation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: warning_continue

implicit none

private
public :: set_increase_factor

    type, public :: Concrete_Adaptation_Parameters
        real(DP) :: increase_factor
        real(DP) :: increase_factor_max
    end type Concrete_Adaptation_Parameters

contains

    subroutine set_increase_factor(classname, current_increase_factor, adaptation_parameters, &
        max_factor_reached)
        character(len=*), intent(in) :: classname
        real(DP), intent(inout) :: current_increase_factor
        type(Concrete_Adaptation_Parameters), intent(in) :: adaptation_parameters
        logical, intent(out) :: max_factor_reached

        current_increase_factor = adaptation_parameters%increase_factor * current_increase_factor
        if (current_increase_factor > adaptation_parameters%increase_factor_max) then
            call warning_continue(classname//": current_increase_factor is too big.")
            current_increase_factor = adaptation_parameters%increase_factor_max
            max_factor_reached = .true.
        else
            max_factor_reached = .false.
        end if
    end subroutine set_increase_factor

end module module_change_adaptation
