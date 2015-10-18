program test_canonical

use, intrinsic :: iso_fortran_env, only: output_unit
use data_constants, only: num_components
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Mixture_Wrapper_Old
use types_mixture_wrapper, only: Mixture_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use types_metropolis_wrapper, only: Metropolis_Wrapper
use procedures_plmc_factory, only: plmc_load, plmc_create, plmc_set, plmc_destroy
use procedures_plmc_visit, only: plmc_visit
use types_observables_wrapper, only: Mixture_Observables_Wrapper
use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use procedures_plmc_propagation, only: plmc_propagator_construct, plmc_propagator_destroy, &
    plmc_propagator_try
use procedures_plmc_write, only: plmc_write
use module_plmc_iterations, only: num_tuning_steps, num_steps, plmc_set_num_steps

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper_Old) :: mixture_old ! to delete
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(num_components)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    type(Mixture_Ewald_Wrapper) :: ewalds
    type(Mixture_Observables_Wrapper) :: observables
    type(Mixture_Observable_Writers_Wrapper) :: observables_writers
    type(Metropolis_Wrapper) :: metropolis

    type(json_file) :: input_data
    integer :: i_step
    logical :: changes_tuned

    call plmc_load(input_data)
    call plmc_set_num_steps(input_data)
    call plmc_create(environment, input_data)
    call plmc_create(mixture, environment, input_data)
    call plmc_create(changes, environment%periodic_box, mixture_old%components, input_data)
    call plmc_create(short_potentials, environment, mixture_old, input_data)
    call plmc_create(ewalds, environment, mixture_old, input_data)
    call plmc_create(observables_writers, environment%walls_potential, mixture_old, changes, input_data)
    call plmc_create(metropolis, environment, changes)
    call input_data%destroy()

    call plmc_set(metropolis, mixture_old%components, short_potentials, ewalds)
    call plmc_propagator_construct(metropolis)
    call plmc_visit(observables, environment%walls_potential, short_potentials, ewalds, mixture_old)
    call plmc_write(-num_tuning_steps, observables_writers, observables)

    if (num_tuning_steps > 0) write(output_unit, *) "Trying to tune changes..."
    do i_step = -num_tuning_steps + 1, 0
        call plmc_propagator_try(metropolis, observables)
        call plmc_set(observables%intras)
        call plmc_set(changes_tuned, i_step, changes, observables%intras)
        call plmc_write(i_step, observables_writers, observables)
        if (changes_tuned) exit
    end do
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call plmc_propagator_try(metropolis, observables)
        call plmc_set(observables%intras)
        call plmc_write(i_step, observables_writers, observables)
    end do
    write(output_unit, *) "Iterations end."

    call plmc_visit(observables, environment%walls_potential, short_potentials, ewalds, mixture_old)
    call plmc_write(i_step-1, observables_writers, observables)

    call plmc_propagator_destroy()
    call plmc_destroy(metropolis)
    call plmc_destroy(observables_writers)
    call plmc_destroy(ewalds)
    call plmc_destroy(short_potentials)
    call plmc_destroy(changes)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program test_canonical
