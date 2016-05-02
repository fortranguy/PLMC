module procedures_plmc_factory

use, intrinsic :: iso_fortran_env, only: error_unit
use data_prefixes, only:mixture_prefix, writers_prefix
use json_module, only: json_file
use procedures_command_arguments, only: create_filename_from_argument
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_physical_model_factory, only: physical_model_create => create, &
    physical_model_destroy => destroy
use procedures_metropolis_algorithms_factory, only: metropolis_algorithms_set
use procedures_markov_chain_generator_factory, only: markov_chain_generator_create => create, &
    markov_chain_generator_destroy => destroy, markov_chain_generator_set => set
use procedures_observables_factory, only: observables_create_generating, &
    observables_destroy_generating
use module_changes_success, only: Concrete_Changes_Success, Concrete_Changes_Counter, &
    Concrete_Switch_Counters, &
    changes_counter_reset, changes_counter_set, switches_counters_reset, switches_counters_set
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_plmc_iterations, only: plmc_set_num_steps
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_readers_wrapper, only: Component_Coordinates_Reader_wrapper, Readers_Wrapper
use procedures_readers_factory, only: readers_create, readers_destroy, &
    readers_set_initial_coordinates
use types_writers_wrapper, only: Writers_Wrapper
use procedures_writers_factory, only: writers_create, writers_destroy

implicit none

private
public :: plmc_create, plmc_destroy, plmc_set

interface plmc_create
    module procedure :: physical_model_create
    module procedure :: markov_chain_generator_create
    module procedure :: observables_create_generating
    module procedure :: create_readers_writers
    module procedure :: create_input_data
end interface plmc_create

interface plmc_destroy
    module procedure :: destroy_input_data
    module procedure :: destroy_readers_writers
    module procedure :: observables_destroy_generating
    module procedure :: markov_chain_generator_destroy
    module procedure :: physical_model_destroy
end interface plmc_destroy

interface plmc_set
    module procedure :: plmc_set_num_steps
    module procedure :: set_mixture_initial_coordinates
    module procedure :: markov_chain_generator_set
    module procedure :: tune_change_components
    module procedure :: set_success_and_reset_counter
end interface plmc_set

contains

    subroutine create_readers_writers(readers, writers, physical_model, changes, input_data)
        type(Readers_Wrapper), intent(out) :: readers
        type(Writers_Wrapper), intent(out) :: writers
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        type(Changes_Wrapper), intent(in) :: changes
        type(json_file), intent(inout) :: input_data

        call readers_create(readers, physical_model%environment%periodic_box, physical_model%&
            mixture%components)
        call writers_create(writers, physical_model%environment, physical_model%short_interactions%&
            wall_pairs, physical_model%mixture%components, changes%components, physical_model%&
            short_interactions%components_pairs, input_data, writers_prefix)
    end subroutine create_readers_writers

    subroutine destroy_readers_writers(readers, writers)
        type(Readers_Wrapper), intent(inout) :: readers
        type(Writers_Wrapper), intent(inout) :: writers

        call writers_destroy(writers)
        call readers_destroy(readers)
    end subroutine destroy_readers_writers

    subroutine create_input_data(input_data, i_argument)
        type(json_file), intent(out) :: input_data
        integer, intent(in) :: i_argument

        character(len=:), allocatable :: data_filename

        call input_data%initialize()
        if (input_data%failed()) call input_data%print_error_message(error_unit)
        call create_filename_from_argument(data_filename, i_argument)
        call input_data%load_file(data_filename)
        if (input_data%failed()) call input_data%print_error_message(error_unit)
    end subroutine create_input_data

    subroutine destroy_input_data(input_data)
        type(json_file), intent(inout) :: input_data

        call input_data%destroy()
    end subroutine destroy_input_data

    subroutine tune_change_components(tuned, i_step, change_components, changes_sucesses)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Component_Wrapper), intent(inout) :: change_components(:)
        type(Concrete_Changes_Success), intent(in) :: changes_sucesses(:)

        logical :: move_tuned(size(change_components)), rotation_tuned(size(change_components))
        integer :: i_component

        do i_component = 1, size(change_components)
            call change_components(i_component)%move_tuner%tune(move_tuned(i_component), i_step, &
                changes_sucesses(i_component)%move)
            call change_components(i_component)%rotation_tuner%tune(rotation_tuned(i_component), &
                i_step, changes_sucesses(i_component)%rotation)
        end do
        tuned = all(move_tuned) .and. all(rotation_tuned)
    end subroutine tune_change_components

    subroutine set_success_and_reset_counter(observables)
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        call changes_counter_set(observables%changes_sucesses, observables%changes_counters)
        call changes_counter_reset(observables%changes_counters)

        call switches_counters_set(observables%switches_successes, observables%switches_counters)
        call switches_counters_reset(observables%switches_counters)
    end subroutine set_success_and_reset_counter

    subroutine set_mixture_initial_coordinates(components_readers, input_data)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_readers(:)
        type(json_file), intent(inout) :: input_data

        call readers_set_initial_coordinates(components_readers, input_data, mixture_prefix)
    end subroutine set_mixture_initial_coordinates

end module procedures_plmc_factory
