program plmc

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use procedures_plmc_propagator, only: plmc_propagator_create, plmc_propagator_destroy, &
    plmc_propagator_try
use types_observables_wrapper, only: Observables_Wrapper
use types_writers_wrapper, only: Writers_Wrapper
use procedures_plmc_factory, only: plmc_create, plmc_set, plmc_destroy
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_visit, only: plmc_visit
use procedures_plmc_write, only: plmc_write
use module_plmc_iterations, only: num_tuning_steps, num_steps, plmc_set_num_steps

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes
    type(Short_Interactions_Wrapper) :: short_interactions
    type(Long_Interactions_Wrapper) :: long_interactions
    type(Metropolis_Algorithms_Wrapper) :: metropolis_algorithms
    type(Observables_Wrapper) :: observables
    type(Writers_Wrapper) :: writers

    type(json_file) :: input_data
    integer :: i_step
    logical :: changes_tuned

    call json_initialize()
    call plmc_create(input_data, 1)
    call plmc_create(environment, input_data)
    call plmc_create(mixture, environment, input_data)
    call plmc_create(short_interactions, environment, mixture, input_data)
    call plmc_create(long_interactions, environment, mixture, input_data)
    call plmc_set_num_steps(input_data)
    call plmc_create(changes, environment%periodic_box, mixture%components, input_data)
    call plmc_create(metropolis_algorithms, environment, changes)
    call plmc_set(metropolis_algorithms, mixture%components, short_interactions, long_interactions)
    call plmc_propagator_create(metropolis_algorithms)
    call plmc_create(observables, mixture%components)
    call plmc_create(writers, mixture%components, short_interactions, long_interactions, changes, &
        input_data)
    call plmc_destroy(input_data)

    call plmc_visit(observables, mixture, short_interactions, long_interactions)
    call plmc_write(-num_tuning_steps, writers, observables)
    if (num_tuning_steps > 0) write(output_unit, *) "Trying to tune changes..."
    do i_step = -num_tuning_steps + 1, 0
        call plmc_propagator_try(observables)
        call plmc_set(observables%changes_sucesses, observables%changes_counters)
        call plmc_set(changes_tuned, i_step, changes%components, observables%changes_sucesses)
        call plmc_write(i_step, writers, observables)
        if (changes_tuned) exit
    end do
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call plmc_propagator_try(observables)
        call plmc_set(observables%changes_sucesses, observables%changes_counters)
        call plmc_write(i_step, writers, observables)
    end do
    write(output_unit, *) "Iterations end."
    call plmc_reset(long_interactions)
    call plmc_visit(observables, mixture, short_interactions, long_interactions)
    call plmc_write(i_step-1, writers, observables)

    call plmc_destroy(writers)
    call plmc_destroy(observables)
    call plmc_propagator_destroy()
    call plmc_destroy(metropolis_algorithms)
    call plmc_destroy(changes)
    call plmc_destroy(long_interactions)
    call plmc_destroy(short_interactions)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program plmc
