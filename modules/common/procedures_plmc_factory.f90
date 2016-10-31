module procedures_plmc_factory

use data_arguments, only: i_generating, i_exploring, num_json_arguments
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_create_output => create_output, json_data_destroy_input => destroy_input, &
    json_data_destroy_output => destroy_output
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: set_nums_particles
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_physical_model_factory, only: physical_model_create_generating => create_generating,&
    physical_model_create_exploring => create_exploring, physical_model_destroy => destroy
use procedures_random_seed_factory, only: random_seed_set => set
use procedures_markov_chain_generator_factory, only: markov_chain_generator_create => create, &
    markov_chain_generator_destroy => destroy
use procedures_markov_chain_explorer_factory, only: markov_chain_explorer_create => create, &
    markov_chain_explorer_destroy => destroy
use procedures_observables_factory, only: observables_create_generating, &
    observables_create_exploring, observables_destroy_generating, observables_destroy_exploring
use module_changes_success, only: Concrete_Changes_Success, reset_counters, set_successes
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_plmc_iterations, only: plmc_set_num_steps, plmc_set_num_snaps
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_readers_wrapper, only: Readers_Wrapper
use procedures_readers_factory, only: readers_create => create, readers_destroy => destroy, &
    readers_set => set
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use procedures_generating_writers_factory, only: generating_writers_create => create, &
    generating_writers_destroy => destroy
use types_exploring_writers_wrapper, only: Exploring_Writers_Wrapper
use procedures_exploring_writers_factory, only: exploring_writers_create => create, &
    exploring_writers_destroy => destroy

implicit none

private
public :: plmc_create, plmc_destroy, plmc_set

interface plmc_create
    module procedure :: physical_model_create_generating, physical_model_create_exploring
    module procedure :: markov_chain_generator_create, markov_chain_explorer_create
    module procedure :: observables_create_generating, observables_create_exploring
    module procedure :: json_data_create_input
    module procedure :: create_generating_data, create_exploring_data
    module procedure :: create_generating_readers_writers, create_exploring_readers_writers
    module procedure :: json_data_create_output
end interface plmc_create

interface plmc_destroy
    module procedure :: json_data_destroy_output
    module procedure :: destroy_generating_readers_writers, destroy_exploring_readers_writers
    module procedure :: destroy_exploring_data
    module procedure :: json_data_destroy_input
    module procedure :: observables_destroy_generating, observables_destroy_exploring
    module procedure :: markov_chain_generator_destroy, markov_chain_explorer_destroy
    module procedure :: physical_model_destroy
end interface plmc_destroy

interface plmc_set
    module procedure :: random_seed_set
    module procedure :: plmc_set_num_steps, plmc_set_num_snaps
    module procedure :: set_initial_observables
    module procedure :: set_coordinates_from_json, set_coordinates_from_snap
    module procedure :: tune_moves
    module procedure :: set_success_and_reset_counter_generating, &
        set_success_and_reset_counter_exploring
end interface plmc_set

contains

    subroutine create_generating_readers_writers(readers, writers, physical_model, changes, &
        generating_data)
        type(Readers_Wrapper), intent(out) :: readers
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes
        type(json_file), intent(inout) :: generating_data

        call readers_create(readers, physical_model%environment, physical_model%mixture%&
            components)
        call generating_writers_create(writers, physical_model%environment, physical_model%&
            short_interactions%wall_pairs, physical_model%mixture%components, physical_model%&
            short_interactions%components_pairs, changes%components, changes%changed_boxes_size, &
            generating_data)
    end subroutine create_generating_readers_writers

    subroutine destroy_generating_readers_writers(readers, writers)
        type(Readers_Wrapper), intent(inout) :: readers
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call generating_writers_destroy(writers)
        call readers_destroy(readers)
    end subroutine destroy_generating_readers_writers

    subroutine create_exploring_readers_writers(readers, writers, physical_model, &
        markov_chain_explorer, visit_energies, generating_data)
        type(Readers_Wrapper), intent(out) :: readers
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Markov_Chain_Explorer_Wrapper), intent(in) :: markov_chain_explorer
        logical, intent(in) :: visit_energies
        type(json_file), intent(inout) :: generating_data

        call readers_create(readers, physical_model%environment, physical_model%mixture%&
            components, markov_chain_explorer%particle_insertion_domains)
        call exploring_writers_create(writers, physical_model%environment, physical_model%&
            short_interactions%wall_pairs, physical_model%mixture%components, physical_model%&
            short_interactions%components_pairs, markov_chain_explorer%volume_change_method, &
            markov_chain_explorer%particle_insertion_method, visit_energies, generating_data)
    end subroutine create_exploring_readers_writers

    subroutine destroy_exploring_readers_writers(readers, writers)
        type(Readers_Wrapper), intent(inout) :: readers
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call readers_destroy(readers)
        call exploring_writers_destroy(writers)
    end subroutine destroy_exploring_readers_writers

    subroutine create_generating_data(generating_data)
        type(json_file), intent(out) :: generating_data

        call json_data_create_input(generating_data, i_generating)
    end subroutine create_generating_data

    subroutine create_exploring_data(generating_data, exploring_data)
        type(json_file), intent(out) :: generating_data, exploring_data

        call json_data_create_input(generating_data, i_generating)
        call json_data_create_input(exploring_data, i_exploring)
    end subroutine create_exploring_data

    subroutine destroy_exploring_data(generating_data, exploring_data)
        type(json_file), intent(inout) :: generating_data, exploring_data

        call json_data_destroy_input(exploring_data)
        call json_data_destroy_input(generating_data)
    end subroutine destroy_exploring_data

    !> @note Beware of inertia
    subroutine tune_moves(tuned, i_step, changes, observables)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Wrapper), intent(inout) :: changes
        type(Generating_Observables_Wrapper), intent(in) :: observables

        logical :: all_boxes_size_tuned,  box_size_tuned
        logical, dimension(size(changes%components, 1), size(changes%components, 2)) :: &
            translation_tuned, rotation_tuned
        integer :: i_box, j_box, i_component

        all_boxes_size_tuned = .true.
        do j_box = 1, size(changes%boxes_size_change_tuner)
            do i_box = 1, size(changes%boxes_size_change_tuner(j_box)%line)
                call changes%boxes_size_change_tuner(j_box)%line(i_box)%tuner%&
                    tune(box_size_tuned, i_step, observables%volumes_change_success(j_box)%&
                        line(i_box))
                all_boxes_size_tuned = all_boxes_size_tuned .and. box_size_tuned
            end do
        end do

        do i_box = 1, size(changes%components, 2)
            do i_component = 1, size(changes%components, 1)
                call changes%components(i_component, i_box)%translation_tuner%&
                    tune(translation_tuned(i_component, i_box), i_step, observables%changes(i_box)%&
                    changes_sucesses(i_component)%translation)
                call changes%components(i_component, i_box)%rotation_tuner%&
                    tune(rotation_tuned(i_component, i_box), i_step, observables%changes(i_box)%&
                    changes_sucesses(i_component)%rotation)
            end do
        end do

        tuned = all_boxes_size_tuned .and. all(translation_tuned) .and. all(rotation_tuned)
    end subroutine tune_moves

    subroutine set_success_and_reset_counter_generating(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        integer :: i_box

        call set_successes(observables%teleportations_successes, observables%&
            teleportations_counters)
        call reset_counters(observables%teleportations_counters)
        call set_successes(observables%volumes_change_success, observables%volumes_change_counter)
        call reset_counters(observables%volumes_change_counter)

        do i_box = 1, size(observables%changes)
            call set_successes(observables%changes(i_box)%changes_sucesses, observables%&
                changes(i_box)%changes_counters)
            call reset_counters(observables%changes(i_box)%changes_counters)
            call set_successes(observables%changes(i_box)%switches_successes, observables%&
                changes(i_box)%switches_counters)
            call reset_counters(observables%changes(i_box)%switches_counters)
            call set_successes(observables%changes(i_box)%transmutations_successes, observables%&
                changes(i_box)%transmutations_counters)
            call reset_counters(observables%changes(i_box)%transmutations_counters)
        end do
    end subroutine set_success_and_reset_counter_generating

    subroutine set_success_and_reset_counter_exploring(observables)
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        integer :: i_box

        do i_box = 1, size(observables%insertion_successes, 2)
            call set_successes(observables%insertion_successes(:, i_box), observables%&
                insertion_counters(:, i_box))
            call reset_counters(observables%insertion_counters(:, i_box))
        end do
    end subroutine set_success_and_reset_counter_exploring

    subroutine set_initial_observables(observables, physical_model)
        type(Generating_Observables_Wrapper), intent(inout) :: observables
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        integer :: i_box

        do i_box = 1, size(observables%accessible_domains_size, 2)
            call set_nums_particles(observables%nums_particles(:, i_box), physical_model%mixture%&
                components(:, i_box))
            observables%accessible_domains_size(:, i_box) = physical_model%environment%&
                accessible_domains(i_box)%get_size()
        end do
    end subroutine set_initial_observables

    subroutine set_coordinates_from_json(readers, generating_data)
        type(Readers_Wrapper), intent(inout) :: readers
        type(json_file), intent(inout) :: generating_data

        call readers_set(readers, generating_data)
    end subroutine set_coordinates_from_json

    subroutine set_coordinates_from_snap(readers, i_snap)
        type(Readers_Wrapper), intent(inout) :: readers
        integer, intent(in) :: i_snap

        call readers_set(readers, num_json_arguments + i_snap)
    end subroutine set_coordinates_from_snap

end module procedures_plmc_factory
