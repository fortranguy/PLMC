module procedures_plmc_write

use data_output_objects, only: random_number_generator_object
use json_module, only: json_core, json_value
use procedures_random_seed_factory, only: random_seed_write => write
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private
public :: plmc_write

interface plmc_write
    module procedure :: write_json_report
    module procedure :: write_generating_observables, write_exploring_observables
end interface plmc_write

contains

    subroutine write_json_report(json, output_data)
        type(json_core), intent(inout) :: json
        type(json_value), intent(inout), pointer :: output_data

        call random_seed_write(json, output_data, random_number_generator_object)
        call json%print(output_data, "report.json")
    end subroutine write_json_report

    subroutine write_generating_observables(writers, observables, num_tuning_steps, num_steps, &
        i_step)
        integer, intent(in) :: num_tuning_steps, num_steps, i_step
        type(Generating_Writers_Wrapper), intent(in) :: writers
        type(Generating_Observables_Wrapper), intent(in) :: observables

        integer :: i_component

        call writers%num_particles%write(i_step, observables%num_particles)
        call writers%field%write(i_step, observables%field_energies)
        call writers%walls%write(i_step, observables%walls_energies)
        if (0 <= i_step) then
            call writers%complete_coordinates%write(i_step)
        end if
        if (-num_tuning_steps < i_step .and. i_step < num_steps) then
            do i_component = 1, size(writers%components_changes)
                call writers%components_changes(i_component)%writer%write(i_step, observables%&
                    changes_sucesses(i_component))
            end do
        end if
        call writers%short_energies%write(i_step, observables%short_energies)
        call writers%dipolar_energies%write(i_step, observables%dipolar_energies)
        call writers%dipolar_mixture_energy%write(i_step, observables%dipolar_mixture_energy)
        call writers%switches%write(i_step, observables%switches_successes)
        call writers%transmutations%write(i_step, observables%transmutations_successes)
    end subroutine write_generating_observables

    subroutine write_exploring_observables(writers, observables, i_snap)
        type(Exploring_Writers_Wrapper), intent(in) :: writers
        type(Exploring_Observables_Wrapper), intent(in) :: observables
        integer, intent(in) :: i_snap

        call writers%beta_pressure_excess%write(i_snap, observables%beta_pressure_excess)
        call writers%field%write(i_snap, observables%field_energies)
        call writers%walls%write(i_snap, observables%walls_energies)
        call writers%inv_pow_activities%write(i_snap, observables%inv_pow_activities)
        call writers%short_energies%write(i_snap, observables%short_energies)
        call writers%dipolar_energies%write(i_snap, observables%dipolar_energies)
        call writers%dipolar_mixture_energy%write(i_snap, observables%dipolar_mixture_energy)
        call writers%insertion_successes%write(i_snap, observables%insertion_successes)
    end subroutine write_exploring_observables

end module procedures_plmc_write
