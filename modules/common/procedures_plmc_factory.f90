module procedures_plmc_factory

use data_input_prefixes, only: environment_prefix, mixture_prefix, random_number_generator_prefix, &
    writers_prefix
use json_module, only: json_file
use procedures_json_data_factory, only: json_data_create_input => create_input, &
    json_data_create_output => create_output, json_data_destroy_input => destroy_input, &
    json_data_destroy_output => destroy_output
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: set_num_particles
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_physical_model_factory, only: physical_model_create_generating => create_generating,&
    physical_model_create_exploring => create_exploring, physical_model_destroy => destroy
use procedures_random_seed_factory, only: random_seed_set => set
use procedures_metropolis_algorithms_factory, only: metropolis_algorithms_set
use procedures_markov_chain_generator_factory, only: markov_chain_generator_create => create, &
    markov_chain_generator_destroy => destroy, markov_chain_generator_set => set
use procedures_markov_chain_explorer_factory, only: markov_chain_explorer_create => create, &
    markov_chain_explorer_destroy => destroy
use procedures_observables_factory, only: observables_create_generating, &
    observables_create_exploring, observables_destroy_generating, observables_destroy_exploring
use module_changes_success, only: Concrete_Changes_Success, &
    reset_counters, set_successes
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_plmc_iterations, only: plmc_set_num_steps, plmc_set_num_snaps
use types_markov_chain_explorer_wrapper, only: Markov_Chain_Explorer_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper
use types_generating_readers_wrapper, only: Generating_Readers_Wrapper
use procedures_generating_readers_factory, only: generating_readers_create => create, &
    generating_readers_destroy => destroy, generating_readers_set => set
use types_exploring_readers_wrapper, only: Exploring_Readers_Wrapper
use procedures_exploring_readers_factory, only: exploring_readers_create => create, &
    exploring_readers_detroy => destroy, exploring_readers_set => set
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
    module procedure :: set_random_seed
    module procedure :: plmc_set_num_steps, set_num_snaps
    module procedure :: set_num_particles
    module procedure :: set_mixture_coordinates, exploring_readers_set
    module procedure :: markov_chain_generator_set
    module procedure :: tune_components_moves
    module procedure :: set_success_and_reset_counter_generating, &
        set_success_and_reset_counter_exploring
end interface plmc_set

contains

    subroutine create_generating_readers_writers(readers, writers, physical_model, changes, &
        generating_data)
        type(Generating_Readers_Wrapper), intent(out) :: readers
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes
        type(json_file), intent(inout) :: generating_data

        call generating_readers_create(readers, physical_model%mixture%components)
        call generating_writers_create(writers, physical_model%environment, physical_model%&
            short_interactions%wall_pairs, physical_model%mixture%components, physical_model%&
            short_interactions%components_pairs, changes%components, generating_data, &
            writers_prefix)
    end subroutine create_generating_readers_writers

    subroutine destroy_generating_readers_writers(readers, writers)
        type(Generating_Readers_Wrapper), intent(inout) :: readers
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call generating_writers_destroy(writers)
        call generating_readers_destroy(readers)
    end subroutine destroy_generating_readers_writers

    subroutine create_exploring_readers_writers(readers, writers, physical_model, num_snaps, &
        markov_chain_explorer, generating_data)
        type(Exploring_Readers_Wrapper), intent(out) :: readers
        type(Exploring_Writers_Wrapper), intent(out) :: writers
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        integer, intent(in) :: num_snaps
        type(Markov_Chain_Explorer_Wrapper), intent(in) :: markov_chain_explorer
        type(json_file), intent(inout) :: generating_data

        integer :: num_offset

        num_offset = 2 ! generating.json + exploring.json, cf. [[create_exploring_data]]
        call exploring_readers_create(readers, physical_model%environment%periodic_box, &
            physical_model%mixture%components, num_snaps, num_offset, generating_data, &
            environment_prefix)
        call exploring_writers_create(writers, physical_model%environment, physical_model%&
            short_interactions%wall_pairs, physical_model%mixture%components, physical_model%&
            short_interactions%components_pairs, markov_chain_explorer%particle_insertion_method)
    end subroutine create_exploring_readers_writers

    subroutine destroy_exploring_readers_writers(readers, writers)
        type(Exploring_Readers_Wrapper), intent(inout) :: readers
        type(Exploring_Writers_Wrapper), intent(inout) :: writers

        call exploring_readers_detroy(readers)
        call exploring_writers_destroy(writers)
    end subroutine destroy_exploring_readers_writers

    subroutine create_generating_data(generating_data)
        type(json_file), intent(out) :: generating_data

        call json_data_create_input(generating_data, 1)
    end subroutine create_generating_data

    subroutine create_exploring_data(generating_data, exploring_data)
        type(json_file), intent(out) :: generating_data, exploring_data

        call json_data_create_input(generating_data, 1)
        call json_data_create_input(exploring_data, 2)
    end subroutine create_exploring_data

    subroutine destroy_exploring_data(generating_data, exploring_data)
        type(json_file), intent(inout) :: generating_data, exploring_data

        call json_data_destroy_input(exploring_data)
        call json_data_destroy_input(generating_data)
    end subroutine destroy_exploring_data

    subroutine set_random_seed(input_data)
        type(json_file), intent(inout) :: input_data

        call random_seed_set(input_data, random_number_generator_prefix)
    end subroutine set_random_seed

    subroutine tune_components_moves(tuned, i_step, changes_components, changes_sucesses)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Component_Wrapper), intent(inout) :: changes_components(:)
        type(Concrete_Changes_Success), intent(in) :: changes_sucesses(:)

        logical :: translation_tuned(size(changes_components)), &
            rotation_tuned(size(changes_components))
        integer :: i_component

        do i_component = 1, size(changes_components)
            call changes_components(i_component)%translation_tuner%&
                tune(translation_tuned(i_component), i_step, changes_sucesses(i_component)%&
                translation)
            call changes_components(i_component)%rotation_tuner%tune(rotation_tuned(i_component), &
                i_step, changes_sucesses(i_component)%rotation)
        end do
        tuned = all(translation_tuned) .and. all(rotation_tuned) ! Beware of inertia.
    end subroutine tune_components_moves

    subroutine set_success_and_reset_counter_generating(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        call set_successes(observables%changes_sucesses, observables%changes_counters)
        call reset_counters(observables%changes_counters)
        call set_successes(observables%switches_successes, observables%switches_counters)
        call reset_counters(observables%switches_counters)
        call set_successes(observables%transmutations_successes, observables%&
            transmutations_counters)
        call reset_counters(observables%transmutations_counters)
    end subroutine set_success_and_reset_counter_generating

    subroutine set_success_and_reset_counter_exploring(observables)
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        call set_successes(observables%insertion_successes, observables%insertion_counters)
        call reset_counters(observables%insertion_counters)
    end subroutine set_success_and_reset_counter_exploring

    subroutine set_mixture_coordinates(components_readers, generating_data)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_readers(:)
        type(json_file), intent(inout) :: generating_data

        call generating_readers_set(components_readers, generating_data, mixture_prefix)
    end subroutine set_mixture_coordinates

    subroutine set_num_snaps(num_snaps, components, generating_data)
        integer, intent(out) :: num_snaps
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: generating_data

        call plmc_set_num_snaps(num_snaps, size(components), generating_data)
    end subroutine set_num_snaps

end module procedures_plmc_factory
