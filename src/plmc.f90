program test_canonical

use, intrinsic :: iso_fortran_env, only: output_unit
use data_constants, only: num_components
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Mixture_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_one_particle_change, only: Abstract_One_Particle_Change
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
use procedures_plmc_tuning, only: plmc_tune, plmc_print_tuning_status

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(num_components)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    type(Mixture_Ewald_Wrapper) :: ewalds
    !class(Abstract_One_Particle_Change), allocatable :: one_particle_move
    !class(Abstract_One_Particle_Rotation), allocatable :: one_particle_rotation
    type(Concrete_Mixture_Observables) :: observables
    type(Mixture_Observable_Writers_Wrapper) :: observables_writers

    type(json_file) :: input_data
    integer :: i_step, num_moves, i_move, num_rotations, i_rotation
    logical :: changes_tuned

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

    !call metropolis_factory_create(one_particle_move, environment, changes)
    !call metropolis_factory_set(one_particle_move, mixture%components, short_potentials, ewalds, &
    !    observables)
    !call metropolis_factory_create(one_particle_rotation, environment, changes)
    !call metropolis_factory_set(one_particle_rotation, mixture%components, ewalds, observables)

    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    num_rotations = mixture%components(1)%orientations%get_num() + &
        mixture%components(2)%orientations%get_num()
    write(output_unit, *)  "Tuning changes..."
    do i_step = 1, num_tuning_steps
        call counters_reset(observables%changes_counters)
        do i_move = 1, num_moves
            !call one_particle_move%try()
        end do
        do i_rotation = 1, num_rotations
            !call one_particle_rotation%try()
        end do
        call sucess_set(observables%changes_success, observables%changes_counters)
        call plmc_write(i_step, observables_writers, observables, in_loop = .true.)
        call plmc_tune(changes_tuned, i_step, changes, observables%changes_success)
        if (changes_tuned) exit
    end do
    call plmc_print_tuning_status(changes_tuned)
    write(output_unit, *) "Iterations start."
    do i_step = 1, num_steps
        call counters_reset(observables%changes_counters)
        do i_move = 1, num_moves
            !call one_particle_move%try()
        end do
        do i_rotation = 1, num_rotations
            !call one_particle_rotation%try()
        end do
        call sucess_set(observables%changes_success, observables%changes_counters)
        call plmc_write(i_step, observables_writers, observables, in_loop = .true.)
    end do
    write(output_unit, *) "Iterations end."

    call plmc_visit(observables, environment%walls_potential, short_potentials, ewalds, mixture)
    call plmc_write(i_step-1, observables_writers, observables, in_loop = .false.)

    !call metropolis_factory_destroy(one_particle_rotation)
    !call metropolis_factory_destroy(one_particle_move)
    call plmc_destroy(observables_writers)
    call plmc_destroy(ewalds)
    call plmc_destroy(short_potentials)
    call plmc_destroy(changes)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program test_canonical
