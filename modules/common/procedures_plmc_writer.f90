module procedures_plmc_writer

use data_output_objects, only: random_number_generator_object
use json_module, only: json_core, json_value
use procedures_random_seed_factory, only: random_seed_write => write
use types_observables_energies, only: Concrete_Observables_Energies
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

        integer :: i_box, j_box, i_component

        do i_box = 1, size(writers%accessible_domains_size)
            call writers%accessible_domains_size(i_box)%write(i_step, observables%&
                accessible_domains_size(:, i_box))
            call writers%nums_particles(i_box)%write(i_step, observables%nums_particles(:, i_box))
            call write_energies(writers%energies(i_box), i_step, observables%energies(i_box))
        end do

        if (-num_tuning_steps < i_step .and. i_step < num_steps) then
            call writers%volumes_change_success%write(i_step, observables%volumes_change_success)
            call writers%volumes_exchange_success%write(i_step, observables%&
                volumes_exchange_success)
            do j_box = 1, size(writers%teleportations_successes, 2)
                do i_box = 1, size(writers%teleportations_successes, 1)
                    call writers%teleportations_successes(i_box, j_box)%writer%&
                        write(i_step, observables%teleportations_successes(:, i_box, j_box))
                    call writers%swaps_successes(i_box, j_box)%writer%&
                        write(i_step, observables%swaps_successes(:, :, i_box, j_box))
                end do
            end do

            do i_box = 1, size(writers%components_changes, 2)
                do i_component = 1, size(writers%components_changes, 1)
                    call writers%components_changes(i_component, i_box)%writer%write(i_step, &
                        observables%changes(i_box)%changes_sucesses(i_component))
                end do
                call writers%switches_successes(i_box)%write(i_step, observables%&
                    changes(i_box)%switches_successes)
                call writers%transmutations_successes(i_box)%write(i_step, observables%&
                    changes(i_box)%transmutations_successes)
            end do
        end if

        if (0 <= i_step) then
            call writers%complete_coordinates%write(i_step)
        end if
    end subroutine write_generating_observables

    subroutine write_exploring_observables(writers, observables, i_snap)
        type(Exploring_Writers_Wrapper), intent(in) :: writers
        type(Exploring_Observables_Wrapper), intent(in) :: observables
        integer, intent(in) :: i_snap

        integer :: i_box

        do i_box = 1, size(writers%maximum_boxes_compression_delta)
            call writers%maximum_boxes_compression_delta(i_box)%write(i_snap, observables%&
                    maximum_boxes_compression_delta(i_box))
            call writers%beta_pressures_excess(i_box)%write(i_snap, observables%&
                beta_pressures_excess(i_box))
            call write_energies(writers%energies(i_box), i_snap, observables%energies(i_box))
            call writers%inv_pow_activities(i_box)%write(i_snap, observables%&
                inv_pow_activities(:, i_box))
            call writers%insertion_successes(i_box)%write(i_snap, observables%&
                insertion_successes(:, i_box))
        end do
    end subroutine write_exploring_observables

    subroutine write_energies(writers, i_step, energies)
        type(Concrete_Energies_Writers), intent(in) :: writers
        integer, intent(in) :: i_step
        type(Concrete_Observables_Energies), intent(in) :: energies

        call writers%field_energies%write(i_step, energies%field_energies)
        call writers%walls_energies%write(i_step, energies%walls_energies)
        call writers%short_energies%write(i_step, energies%short_energies)
        call writers%dipolar_energies%write(i_step, energies%dipolar_energies)
        call writers%dipolar_shared_energy%write(i_step, energies%dipolar_shared_energy)
    end subroutine write_energies

end module procedures_plmc_writer
