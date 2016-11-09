module procedures_markov_chain_generator_factory

use json_module, only: json_file
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_changes_factory, only: changes_create => create, changes_destroy => destroy
use procedures_generating_algorithms_factory, only: generating_algorithms_create => create, &
    generating_algorithms_destroy => destroy
use procedures_plmc_propagator_factory, only: plmc_propagator_create => create, &
    plmc_propagator_destroy => destroy
use types_markov_chain_generator_wrapper, only: Markov_Chain_Generator_Wrapper

implicit none

private
public :: create, destroy

contains

    subroutine create(markov_chain_generator, physical_model, num_tuning_steps, generating_data)
        type(Markov_Chain_Generator_Wrapper), intent(out) :: markov_chain_generator
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data

        call changes_create(markov_chain_generator%changes, physical_model%environment, &
            physical_model%mixture%components, num_tuning_steps, generating_data)
        call generating_algorithms_create(markov_chain_generator%generating_algorithms, &
            physical_model, markov_chain_generator%changes)
        call plmc_propagator_create(markov_chain_generator%plmc_propagator, physical_model%mixture%&
            components, markov_chain_generator%generating_algorithms, generating_data)
    end subroutine create

    subroutine destroy(markov_chain_generator)
        type(Markov_Chain_Generator_Wrapper), intent(inout) :: markov_chain_generator

        call plmc_propagator_destroy(markov_chain_generator%plmc_propagator)
        call generating_algorithms_destroy(markov_chain_generator%generating_algorithms)
        call changes_destroy(markov_chain_generator%changes)
    end subroutine destroy

end module procedures_markov_chain_generator_factory
