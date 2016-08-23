program plmc_generate

use, intrinsic :: iso_fortran_env, only: output_unit
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_markov_chain_generator_wrapper, only: Markov_Chain_Generator_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_generating_io, only: Generating_IO_Wrapper
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_visit, only: plmc_visit
use procedures_plmc_write, only: plmc_write
use procedures_plmc_help, only: plmc_catch_generating_help

implicit none

    type(Physical_Model_Wrapper) :: physical_model
    type(Markov_Chain_Generator_Wrapper) :: markov_chain_generator
    type(Generating_Observables_Wrapper) :: observables
    type(Generating_IO_Wrapper) :: io

    integer :: num_tuning_steps, num_steps, i_step
    logical :: changes_tuned

    call plmc_catch_generating_help()

    call plmc_create(io%generating_data)
    call plmc_create(physical_model, io%generating_data)
    call plmc_set(io%generating_data)
    call plmc_set(num_tuning_steps, num_steps, io%generating_data)
    call plmc_create(markov_chain_generator, physical_model, num_tuning_steps, io%generating_data)
    call plmc_create(observables, physical_model%mixture%components)
    call plmc_create(io%readers, io%writers, physical_model, markov_chain_generator%changes, &
        io%generating_data)
    call plmc_set(io%readers, io%generating_data)
    call plmc_set(markov_chain_generator)
    call plmc_destroy(io%generating_data)
    call plmc_create(io%json, io%report_data)
    call plmc_write(io%json, io%report_data)
    call plmc_destroy(io%json, io%report_data)

    call plmc_reset(physical_model)
    call plmc_set(observables, physical_model) !in exploring too?
    call plmc_visit(observables%energies, physical_model)
    call plmc_write(io%writers, observables, num_tuning_steps, num_steps, -num_tuning_steps)
    if (num_tuning_steps > 0) write(output_unit, *) "Trying to tune changes..."
    do i_step = -num_tuning_steps + 1, 0
        call markov_chain_generator%plmc_propagator%try(observables)
        call plmc_set(observables)
        call plmc_set(changes_tuned, i_step, markov_chain_generator%changes%components, &
            observables%changes_sucesses)
        call plmc_write(io%writers, observables, num_tuning_steps, num_steps, i_step)
        if (changes_tuned) exit
    end do
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call markov_chain_generator%plmc_propagator%try(observables)
        call plmc_set(observables)
        call plmc_write(io%writers, observables, num_tuning_steps, num_steps, i_step)
    end do
    write(output_unit, *) "Iterations end."
    call plmc_reset(physical_model)
    call plmc_visit(observables%energies, physical_model)
    call plmc_write(io%writers, observables, num_tuning_steps, num_steps, num_steps)

    call plmc_destroy(io%readers, io%writers)
    call plmc_destroy(observables)
    call plmc_destroy(markov_chain_generator)
    call plmc_destroy(physical_model)

end program plmc_generate
