module procedures_plmc_tune

use types_changes_wrapper, only: Mixture_Changes_Wrapper
use types_mixture_observables, only: Concrete_Mixture_Observables

implicit none

private
public :: plmc_tune

contains

    subroutine plmc_tune(tuned, i_step, changes, observables)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Mixture_Changes_Wrapper), intent(inout) :: changes
        type(Concrete_Mixture_Observables), intent(in) :: observables

        tuned = .false.
        !call changes%move_tuner%tune(move_tune, i_step, observables%changes_success%ratios%move)
    end subroutine plmc_tune

end module procedures_plmc_tune
