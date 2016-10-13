module procedures_markov_chain_generator_factory

use json_module, only: json_file
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_changes_factory, only: changes_create, changes_destroy
use procedures_generating_algorithms_factory, only: generating_algorithms_create, &
    generating_algorithms_destroy, generating_algorithms_set
use procedures_plmc_propagator_factory, only: plmc_propagator_create => create, &
    plmc_propagator_destroy => destroy
use types_markov_chain_generator_wrapper, only: Markov_Chain_Generator_Wrapper

implicit none

private
public :: create, destroy, set

contains

    subroutine create(markov_chain_generator, physical_model, num_tuning_steps, generating_data)
        type(Markov_Chain_Generator_Wrapper), intent(out) :: markov_chain_generator
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data

        call changes_create(markov_chain_generator%changes, physical_model%environment, &
            physical_model%mixture%gemc_components, num_tuning_steps, generating_data)
        call generating_algorithms_create(markov_chain_generator%generating_algorithms, &
            physical_model, markov_chain_generator%changes)
        call plmc_propagator_create(markov_chain_generator%plmc_propagator, markov_chain_generator%&
            generating_algorithms)
    end subroutine create

    subroutine destroy(markov_chain_generator)
        type(Markov_Chain_Generator_Wrapper), intent(inout) :: markov_chain_generator

        call plmc_propagator_destroy(markov_chain_generator%plmc_propagator)
        call generating_algorithms_destroy(markov_chain_generator%generating_algorithms)
        call changes_destroy(markov_chain_generator%changes)
    end subroutine destroy

    subroutine set(markov_chain_generator)
        type(Markov_Chain_Generator_Wrapper), intent(inout) :: markov_chain_generator

        call generating_algorithms_set(markov_chain_generator%generating_algorithms)
        call markov_chain_generator%plmc_propagator%set_selector()
    end subroutine set

end module procedures_markov_chain_generator_factory
