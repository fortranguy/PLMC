module procedures_plmc_write

use data_output_objects, only: random_number_generator_object
use json_module, only: json_core, json_value
use procedures_random_seed_factory, only: random_seed_write => write
use types_observables_energies, only: Concrete_Energies
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_energies_writers, only: Concrete_Energies_Writers
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper

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

        call writers%accessible_domain_size%write(i_step, observables%accessible_domain_size)
        call writers%volume_change_success%write(i_step, observables%volume_change_success)
        call writers%nums_particles%write(i_step, observables%nums_particles)
        if (0 <= i_step) then
            call writers%complete_coordinates%write(i_step)
        end if
        call write_energies(writers%energies, observables%energies, i_step)
        if (-num_tuning_steps < i_step .and. i_step < num_steps) then
            do i_component = 1, size(writers%components_changes)
                call writers%components_changes(i_component)%writer%write(i_step, observables%&
                    changes_sucesses(i_component))
            end do
        end if
        call writers%switches_successes%write(i_step, observables%switches_successes)
        call writers%transmutations_successes%write(i_step, observables%transmutations_successes)
    end subroutine write_generating_observables

    subroutine write_exploring_observables(writers, observables, i_snap)
        type(Exploring_Writers_Wrapper), intent(in) :: writers
        type(Exploring_Observables_Wrapper), intent(in) :: observables
        integer, intent(in) :: i_snap

        call writers%maximum_box_compression_delta%write(i_snap, observables%&
            maximum_box_compression_delta)
        call writers%beta_pressure_excess%write(i_snap, observables%beta_pressure_excess)
        call write_energies(writers%energies, observables%energies, i_snap)
        call writers%inv_pow_activities%write(i_snap, observables%inv_pow_activities)
        call writers%insertion_successes%write(i_snap, observables%insertion_successes)
    end subroutine write_exploring_observables

    subroutine write_energies(writers, energies, i_step)
        type(Concrete_Energies_Writers), intent(in) :: writers
        type(Concrete_Energies), intent(in) :: energies
        integer, intent(in) :: i_step

        call writers%field_energies%write(i_step, energies%field_energies)
        call writers%walls_energies%write(i_step, energies%walls_energies)
        call writers%short_energies%write(i_step, energies%short_energies)
        call writers%dipolar_energies%write(i_step, energies%dipolar_energies)
        call writers%dipolar_shared_energy%write(i_step, energies%dipolar_shared_energy)
    end subroutine write_energies

end module procedures_plmc_write
