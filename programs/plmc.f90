program plmc

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_markov_chain_generator_wrapper, only: Markov_Chain_Generator_Wrapper
use types_input_output_wrapper, only: Input_Output_Wrapper

use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator
use types_observables_wrapper, only: Observables_Wrapper
use types_readers_wrapper, only: Readers_Wrapper
use types_writers_wrapper, only: Writers_Wrapper
use module_plmc_iterations, only: num_tuning_steps, num_steps
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_visit, only: plmc_visit
use procedures_plmc_write, only: plmc_write
use procedures_plmc_help, only: plmc_catch_help

implicit none

    type(Physical_Model_Wrapper) :: physics
    type(Markov_Chain_Generator_Wrapper) :: markov_chain
    type(Observables_Wrapper) :: observables
    type(Input_Output_Wrapper) :: io

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Short_Interactions_Wrapper) :: short_interactions
    type(Dipolar_Interactions_Wrapper) :: dipolar_interactions
    type(Changes_Wrapper) :: changes
    type(Metropolis_Algorithms_Wrapper) :: metropolis_algorithms
    class(Abstract_PLMC_Propagator), allocatable :: plmc_propagator

    type(Readers_Wrapper) :: readers
    type(Writers_Wrapper) :: writers

    type(json_file) :: input_data
    integer :: i_step
    logical :: changes_tuned

    call plmc_catch_help()
    call plmc_create(input_data, 1)
    call plmc_create(environment, mixture, short_interactions, dipolar_interactions, changes, &
        metropolis_algorithms, plmc_propagator, observables, readers, writers, input_data)
    call plmc_set(readers%components, input_data)
    call plmc_set(metropolis_algorithms)
    call plmc_propagator%set_selector()
    call plmc_destroy(input_data)

    call plmc_reset(mixture%total_moment, short_interactions, dipolar_interactions)
    call plmc_visit(observables, environment, mixture, short_interactions, dipolar_interactions)
    call plmc_write(-num_tuning_steps, writers, observables)
    if (num_tuning_steps > 0) write(output_unit, *) "Trying to tune changes..."
    do i_step = -num_tuning_steps + 1, 0
        call plmc_propagator%try(observables)
        call plmc_set(observables)
        call plmc_set(changes_tuned, i_step, changes%components, observables%changes_sucesses)
        call plmc_write(i_step, writers, observables)
        if (changes_tuned) exit
    end do
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call plmc_propagator%try(observables)
        call plmc_set(observables)
        call plmc_write(i_step, writers, observables)
    end do
    write(output_unit, *) "Iterations end."
    call plmc_reset(mixture%total_moment, short_interactions, dipolar_interactions)
    call plmc_visit(observables, environment, mixture, short_interactions, dipolar_interactions)
    call plmc_write(i_step-1, writers, observables)

    call plmc_destroy(environment, mixture, short_interactions, dipolar_interactions, changes, &
        metropolis_algorithms, plmc_propagator, observables, readers, writers)

end program plmc
