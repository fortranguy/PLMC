!> @note After thermalisation, plmc_propagator%reset() was removed because it may delay
!> the search of the equilibrium (e.g. GEMC: DHS+HS at low temperature).
program plmc_generate

use, intrinsic :: iso_fortran_env, only: output_unit
use data_output_objects, only: generating_report_filename
use data_arguments, only: i_generating
use json_module, only: json_core
use procedures_json_data_factory, only: json_data_create => create, json_data_destroy => destroy
use procedures_json_reports_factory, only: json_reports_create => create, json_reports_destroy => &
    destroy
use procedures_random_seed_factory, only: random_seed_add_to_report => add_to_report
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use procedures_generating_algorithms_factory, only: generating_algorithms_add_to_report => &
    add_to_report
use types_markov_chain_generator_wrapper, only: Markov_Chain_Generator_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_generating_io, only: Generating_IO_Wrapper
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_resetter, only: plmc_reset
use procedures_plmc_visitor, only: plmc_visit
use procedures_plmc_writer, only: plmc_write
use procedures_plmc_help, only: plmc_catch_generating_help

implicit none

    type(Physical_Model_Wrapper) :: physical_model
    type(Markov_Chain_Generator_Wrapper) :: markov_chain_generator
    type(Generating_Observables_Wrapper) :: observables
    type(json_core) :: json
    type(Generating_IO_Wrapper) :: io

    integer :: num_tuning_steps, num_steps, i_step
    logical :: changes_tuned

    call plmc_catch_generating_help()

    call json_data_create(io%data, i_generating)
    call plmc_create(physical_model, io%data)
    call plmc_set(io%data)
    call plmc_set(num_tuning_steps, num_steps, io%data)
    call plmc_create(markov_chain_generator, physical_model, num_tuning_steps, io%data)
    call plmc_create(observables, physical_model%mixture%components)
    call plmc_create(io%readers, io%writers, physical_model, markov_chain_generator%changes, io%&
        data)
    call plmc_set(io%readers, io%data)
    call json_data_destroy(io%data)
    call json%initialize()
    call json%create_object(io%report%root, "")
    call json_reports_create(json, io%report)
    call random_seed_add_to_report(json, io%report%random_seed, "initial seed")
    call json%print(io%report%root, generating_report_filename)

    call plmc_reset(physical_model, skip_dipolar_interactions=.false.)
    call markov_chain_generator%plmc_propagator%reset()
    call generating_algorithms_add_to_report(json, io%report%algorithms_weight, &
        markov_chain_generator%generating_algorithms)
    call json%print(io%report%root, generating_report_filename)
    call plmc_set(observables, physical_model) !in exploring too?
    call plmc_visit(observables%energies, physical_model, use_cells=.false.)
    call plmc_write(io%writers, observables, num_tuning_steps, num_steps, -num_tuning_steps)
    if (num_tuning_steps > 0) write(output_unit, *) "Trying to tune changes..."
    do i_step = -num_tuning_steps + 1, 0
        call markov_chain_generator%plmc_propagator%try(observables)
        call plmc_set(changes_tuned, i_step, markov_chain_generator%changes, observables)
        call plmc_set(observables)
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
    call plmc_reset(physical_model, skip_dipolar_interactions=.false.)
    call plmc_visit(observables%energies, physical_model, use_cells=.false.)
    call plmc_write(io%writers, observables, num_tuning_steps, num_steps, num_steps)
    call random_seed_add_to_report(json, io%report%random_seed, "final seed")
    call json%print(io%report%root, generating_report_filename)

    call json_reports_destroy(json, io%report)
    call plmc_destroy(io%readers, io%writers)
    call plmc_destroy(observables)
    call plmc_destroy(markov_chain_generator)
    call plmc_destroy(physical_model)

end program plmc_generate
