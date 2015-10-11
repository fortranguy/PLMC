module procedures_plmc_tuning

use, intrinsic :: iso_fortran_env, only: output_unit
use data_constants, only: num_components
use module_plmc_iterations, only: num_tuning_steps
use types_changes_wrapper, only: Changes_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: plmc_tune

contains

    subroutine plmc_tune(tuned, i_step, changes, observables_intras)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Wrapper), intent(inout) :: changes(:)
        type(Observables_Wrapper), intent(in) :: observables_intras(num_components)

        logical :: move_tuned(num_components), rotation_tuned(num_components)
        integer :: i_component

        do i_component = 1, num_components
            call changes(i_component)%move_tuner%tune(move_tuned(i_component), i_step, &
                observables_intras(i_component)%changes_success%move)
            call changes(i_component)%rotation_tuner%tune(rotation_tuned(i_component), i_step, &
                observables_intras(i_component)%changes_success%rotation)
        end do
        tuned = all(move_tuned) .and. all(rotation_tuned)
    end subroutine plmc_tune

end module procedures_plmc_tuning
