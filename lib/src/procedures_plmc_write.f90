module procedures_plmc_write

use module_plmc_iterations, only: num_tuning_steps, num_steps
use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use types_observables_wrapper, only: Mixture_Observables_Wrapper

implicit none

private
public :: plmc_write

interface plmc_write
    module procedure :: write_observables
end interface plmc_write

contains

    subroutine write_observables(i_step, observables_writers, observables)
        integer, intent(in) :: i_step
        type(Mixture_Observable_Writers_Wrapper), intent(in) :: observables_writers
        type(Mixture_Observables_Wrapper), intent(in) :: observables

        if (i_step > 0) then
            call observables_writers%intras(1)%coordinates%write(i_step)
            call observables_writers%intras(2)%coordinates%write(i_step)
        end if
        call observables_writers%intras(1)%energy%write(i_step, &
            observables%intras(1)%particles_energy)
        call observables_writers%intras(2)%energy%write(i_step, &
            observables%intras(2)%particles_energy)
        call observables_writers%inter_energy%write(i_step, observables%inter_energy)
        if (-num_tuning_steps < i_step .and. i_step < num_steps) then
            call observables_writers%intras(1)%changes%write(i_step, &
                observables%intras(1)%changes_success)
            call observables_writers%intras(2)%changes%write(i_step, &
                observables%intras(2)%changes_success)
        end if
    end subroutine write_observables

end module procedures_plmc_write
