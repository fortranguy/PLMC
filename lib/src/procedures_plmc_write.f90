module procedures_plmc_write

use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use types_mixture_observables, only: Concrete_Mixture_Observables

implicit none

private
public :: plmc_write

interface plmc_write
    module procedure :: write_observables
end interface plmc_write

contains

    subroutine write_observables(i_step, observables_writers, observables, in_loop)
        integer, intent(in) :: i_step
        type(Mixture_Observable_Writers_Wrapper), intent(in) :: observables_writers
        type(Concrete_Mixture_Observables), intent(in) :: observables
        logical, intent(in) :: in_loop

        call observables_writers%intras(1)%coordinates%write(i_step)
        call observables_writers%intras(2)%coordinates%write(i_step)
        call observables_writers%intras(1)%energy%write(i_step, observables%particles_energies(1))
        call observables_writers%intras(2)%energy%write(i_step, observables%particles_energies(2))
        call observables_writers%inter_energy%write(i_step, observables%inter_energy)
        if (in_loop) then
            call observables_writers%intras(1)%changes%write(i_step, &
                observables%changes_success%ratios(1))
            call observables_writers%intras(2)%changes%write(i_step, &
                observables%changes_success%ratios(2))
        end if
    end subroutine write_observables

end module procedures_plmc_write
