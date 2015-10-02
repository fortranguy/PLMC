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
use module_particles_energy, only: Concrete_Particles_Energy, &
    particle_energys_write_legend => Concrete_Particles_Energy_write_legend, &
    particles_energy_write => Concrete_Particles_Energy_write
use procedures_visits, only: visit
use procedures_plmc_factory, only: plmc_load, plmc_create, plmc_destroy
use types_mixture_observables, only: Concrete_Mixture_Observables

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(2)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Mixture_Observables) :: observables

    type(json_file) :: input_data
    character(len=:), allocatable :: data_field
    logical :: data_found
    integer :: energy_units(2), inter_energy_unit
    integer :: num_steps, i_step, num_moves, i_move
    integer :: move_units(2)

    call plmc_load(input_data)
    call plmc_create(environment, input_data)
    call plmc_create(mixture, input_data, environment)
    call plmc_create(changes, input_data, environment%periodic_box, mixture%components)
    call plmc_create(short_potentials, input_data, environment, mixture)

    data_field = "Monte Carlo.number of steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)
    call check_positive("test_canonical", "num_steps", num_steps)
    deallocate(data_field)

    call input_data%destroy()

    call visit(observables%particles_energies, observables%inter_energy, short_potentials%intras, &
        short_potentials%inter_micro%pair, mixture%components)
    call visit(observables%particles_energies, environment%walls_potential, mixture%components, &
        short_potentials%walls)

    i_step = 0
    open(newunit=energy_units(1), recl=4096, file="component_1_energy.out", action="write")
    call particle_energys_write_legend(energy_units(1))
    call particles_energy_write(energy_units(1), i_step, observables%particles_energies(1))
    open(newunit=energy_units(2), recl=4096, file="component_2_energy.out", action="write")
    call particle_energys_write_legend(energy_units(2))
    call particles_energy_write(energy_units(2), i_step, observables%particles_energies(2))
    open(newunit=inter_energy_unit, recl=4096, file="inter_12_energy.out", action="write")
    write(inter_energy_unit, *) "# i_step    inter"
    write(inter_energy_unit, *) i_step, observables%inter_energy

    call metropolis_factory_create(one_particle_move, environment, changes)
    call metropolis_factory_set(one_particle_move, mixture%components)
    call metropolis_factory_set(one_particle_move, short_potentials%intras, &
        short_potentials%inters, short_potentials%walls)
    call metropolis_factory_set(one_particle_move, observables)

    open(newunit=move_units(1), recl=4096, file="component_1_move.out", action="write")
    write(move_units(1), *) "# num_hits    num_success"
    open(newunit=move_units(2), recl=4096, file="component_2_move.out", action="write")
    write(move_units(2), *) "# num_hits    num_success"

    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    do i_step = 1, num_steps
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call particles_energy_write(energy_units(1), i_step, observables%particles_energies(1))
        call particles_energy_write(energy_units(2), i_step, observables%particles_energies(2))
        write(inter_energy_unit, *) i_step, observables%inter_energy
        write(move_units(1), *) i_step, observables%move_counters(1)%num_hits, observables%move_counters(1)%num_success
        write(move_units(2), *) i_step, observables%move_counters(2)%num_hits, observables%move_counters(2)%num_success
        observables%move_counters(1)%num_hits = 0; observables%move_counters(1)%num_success = 0
        observables%move_counters(2)%num_hits = 0; observables%move_counters(2)%num_success = 0
    end do

    call visit(observables%particles_energies, observables%inter_energy, short_potentials%intras, &
        short_potentials%inter_micro%pair, mixture%components)
    call visit(observables%particles_energies, environment%walls_potential, mixture%components, &
        short_potentials%walls)
    call particles_energy_write(energy_units(1), i_step-1, observables%particles_energies(1))
    call particles_energy_write(energy_units(2), i_step-1, observables%particles_energies(2))
    write(inter_energy_unit, *) i_step-1, observables%inter_energy

    close(move_units(1))
    close(move_units(2))
    call metropolis_factory_destroy(one_particle_move)
    close(inter_energy_unit)
    close(energy_units(2))
    close(energy_units(1))
    call plmc_destroy(short_potentials)
    call plmc_destroy(changes)
    call plmc_destroy(mixture)
    call plmc_destroy(environment)

end program test_canonical
