module procedures_plmc_write

use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private
public :: plmc_write

interface plmc_write
    module procedure :: write_generating_observables, write_exploring_observables
end interface plmc_write

contains

    subroutine write_generating_observables(writers, observables, num_tuning_steps, num_steps, &
        i_step)
        integer, intent(in) :: num_tuning_steps, num_steps, i_step
        type(Generating_Writers_Wrapper), intent(in) :: writers
        type(Generating_Observables_Wrapper), intent(in) :: observables

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
    end subroutine write_generating_observables

    subroutine write_exploring_observables(writers, observables, i_snap)
        type(Exploring_Writers_Wrapper), intent(in) :: writers
        type(Exploring_Observables_Wrapper), intent(in) :: observables
        integer, intent(in) :: i_snap

        call writers%widom_successes%write(i_snap, observables%widom_successes)
        call writers%inv_pow_activities%write(i_snap, observables%inv_pow_activities)
    end subroutine write_exploring_observables

end module procedures_plmc_write
