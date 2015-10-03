program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_checks, only: check_positive
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Mixture_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use class_one_particle_move, only: Abstract_One_Particle_Move
use procedures_metropolis_factory, only: metropolis_factory_create, metropolis_factory_set, &
    metropolis_factory_destroy
use module_particles_energy, only: Concrete_Particles_Energy
use procedures_plmc_factory, only: plmc_load, plmc_create, plmc_destroy
use procedures_plmc_visit, only: plmc_visit
use module_changes_success, only: counters_reset => Concrete_Mixture_Changes_Counters_reset, &
    sucess_set => Concrete_Mixture_Changes_Success_set
use types_mixture_observables, only: Concrete_Mixture_Observables
use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use procedures_plmc_write, only: plmc_write

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(2)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Mixture_Observables) :: observables
    type(Mixture_Observable_Writers_Wrapper) :: observables_writers

    type(json_file) :: input_data
    character(len=:), allocatable :: data_field
    logical :: data_found
    integer :: num_steps, i_step, num_moves, i_move

    call plmc_load(input_data)
    call plmc_create(environment, input_data)
    call plmc_create(mixture, input_data, environment)
    call plmc_create(changes, input_data, environment%periodic_box, mixture%components)
    call plmc_create(short_potentials, input_data, environment, mixture)
    call plmc_create(observables_writers, mixture, changes)

    data_field = "Monte Carlo.number of steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)
    call check_positive("test_canonical", "num_steps", num_steps)
    deallocate(data_field)

    call input_data%destroy()

    call plmc_visit(observables, environment%walls_potential, short_potentials, mixture)
    call plmc_write(0, observables_writers, observables, in_loop = .false.)

    call metropolis_factory_create(one_particle_move, environment, changes)
    call metropolis_factory_set(one_particle_move, mixture%components)
    call metropolis_factory_set(one_particle_move, short_potentials%intras, short_potentials%inters)
    call metropolis_factory_set(one_particle_move, observables)

    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    do i_step = 1, num_steps
        call counters_reset(observables%changes_counters)
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call sucess_set(observables%changes_success, observables%changes_counters)
        call plmc_write(i_step, observables_writers, observables, in_loop = .true.)
    end do

    call plmc_visit(observables, environment%walls_potential, short_potentials, mixture)
    call plmc_write(i_step-1, observables_writers, observables, in_loop = .false.)

    call metropolis_factory_destroy(one_particle_move)
    call plmc_destroy(observables_writers)
    call plmc_destroy(short_potentials)
    call plmc_destroy(changes)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program test_canonical
