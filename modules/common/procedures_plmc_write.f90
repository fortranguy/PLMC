module procedures_plmc_write

use types_writers_wrapper, only: Writers_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: plmc_write

interface plmc_write
    module procedure :: write_observables
end interface plmc_write

contains

    subroutine write_observables(num_tuning_steps, num_steps, i_step, writers, observables)
        integer, intent(in) :: num_tuning_steps, num_steps, i_step
        type(Writers_Wrapper), intent(in) :: writers
        type(Observables_Wrapper), intent(in) :: observables

        integer :: i_component

        if (i_step > 0) then
            do i_component = 1, size(writers%components)
                call writers%components(i_component)%coordinates%write(i_step)
            end do
        end if
        if (-num_tuning_steps < i_step .and. i_step < num_steps) then
            do i_component = 1, size(writers%components)
                call writers%components(i_component)%changes%write(i_step, observables%&
                    changes_sucesses(i_component))
            end do
        end if
        call writers%field%write(i_step, observables%field_energies)
        call writers%walls%write(i_step, observables%walls_energies)
        call writers%switches%write(i_step, observables%switches_successes)
        call writers%short_energies%write(i_step, observables%short_energies)
        call writers%dipolar_energies%write(i_step, observables%dipolar_energies)
        call writers%dipolar_mixture_energy%write(i_step, observables%dipolar_mixture_energy)
    end subroutine write_observables

end module procedures_plmc_write
