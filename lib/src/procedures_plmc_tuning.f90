module procedures_plmc_tuning

use, intrinsic :: iso_fortran_env, only: output_unit
use data_constants, only: num_components
use module_plmc_iterations, only: num_tuning_steps
use types_changes_wrapper, only: Changes_Wrapper
use module_changes_success, only: Concrete_Mixture_Changes_Success

implicit none

private
public :: plmc_tune

contains

    subroutine plmc_tune(tuning_is_over, i_step, changes, changes_success)
        logical, intent(out) :: tuning_is_over
        integer, intent(in) :: i_step
        type(Changes_Wrapper), intent(inout) :: changes(:)
        type(Concrete_Mixture_Changes_Success), intent(in) :: changes_success

        logical :: tuned

        tuning_is_over = num_tuning_steps < i_step
        if (.not. tuning_is_over) then
            call tune_changes(tuned, i_step, changes, changes_success)
            if (tuned) then
                write(output_unit, *) "Changes may be tuned."
                tuning_is_over = .true.
                return
            end if
        else
            write(output_unit, *) "Changes may not be tuned."
        end if
    end subroutine plmc_tune

    subroutine tune_changes(tuned, i_step, changes, changes_success)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Wrapper), intent(inout) :: changes(:)
        type(Concrete_Mixture_Changes_Success), intent(in) :: changes_success

        logical :: move_tuned(num_components), rotation_tuned(num_components)
        integer :: i_component

        do i_component = 1, num_components
            call changes(i_component)%move_tuner%tune(move_tuned(i_component), i_step, &
                changes_success%ratios(i_component)%move)
            call changes(i_component)%rotation_tuner%tune(rotation_tuned(i_component), i_step, &
                changes_success%ratios(i_component)%rotation)
        end do
        tuned = all(move_tuned) .and. all(rotation_tuned)
    end subroutine tune_changes

end module procedures_plmc_tuning
