program test_canonical

use, intrinsic :: iso_fortran_env, only: output_unit
use data_constants, only: num_components
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Mixture_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_one_particle_move, only: Abstract_One_Particle_Move
use procedures_metropolis_factory, only: metropolis_factory_create, metropolis_factory_set, &
    metropolis_factory_destroy
use procedures_plmc_factory, only: plmc_load, plmc_create, plmc_destroy
use procedures_plmc_visit, only: plmc_visit
use module_changes_success, only: counters_reset => Concrete_Mixture_Changes_Counters_reset, &
    sucess_set => Concrete_Mixture_Changes_Success_set
use types_mixture_observables, only: Concrete_Mixture_Observables
use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use procedures_plmc_write, only: plmc_write
use module_plmc_iterations, only: num_tuning_steps, num_steps, plmc_set_num_steps
use procedures_plmc_tuning, only: plmc_tune

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(num_components)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    type(Mixture_Ewald_Wrapper) :: ewalds
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Mixture_Observables) :: observables
    type(Mixture_Observable_Writers_Wrapper) :: observables_writers

    type(json_file) :: input_data
    integer :: i_step, num_moves, i_move
    logical :: tuning_is_over

    call plmc_load(input_data)
    call plmc_set_num_steps(input_data)
    call plmc_create(environment, input_data)
    call plmc_create(mixture, environment, input_data)
    call plmc_create(changes, environment%periodic_box, mixture%components, input_data)
    call plmc_create(short_potentials, environment, mixture, input_data)
    call plmc_create(ewalds, environment, mixture, input_data)
    call plmc_create(observables_writers, environment%walls_potential, mixture, changes, input_data)
    call input_data%destroy()

    call plmc_visit(observables, environment%walls_potential, short_potentials, ewalds, mixture)
    call plmc_write(0, observables_writers, observables, in_loop = .false.)

    call metropolis_factory_create(one_particle_move, environment, changes)
    call metropolis_factory_set(one_particle_move, mixture%components, short_potentials, ewalds, &
        observables)

    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    write(output_unit, *)  "Tuning changes..."
    i_step = 1; tuning_is_over = .false.
    do while (.not. tuning_is_over)
        call counters_reset(observables%changes_counters)
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call sucess_set(observables%changes_success, observables%changes_counters)
        call plmc_tune(tuning_is_over, i_step, changes, observables%changes_success)
        call plmc_write(i_step, observables_writers, observables, in_loop = .true.)
        i_step = i_step + 1
    end do
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call counters_reset(observables%changes_counters)
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call sucess_set(observables%changes_success, observables%changes_counters)
        call plmc_write(i_step, observables_writers, observables, in_loop = .true.)
    end do
    write(output_unit, *) "Iterations end."

    call plmc_visit(observables, environment%walls_potential, short_potentials, ewalds, mixture)
    call plmc_write(i_step-1, observables_writers, observables, in_loop = .false.)

    call metropolis_factory_destroy(one_particle_move)
    call plmc_destroy(observables_writers)
    call plmc_destroy(ewalds)
    call plmc_destroy(short_potentials)
    call plmc_destroy(changes)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program test_canonical
