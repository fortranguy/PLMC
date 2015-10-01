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
use types_change_counter, only: Concrete_Change_Counter
use module_particle_energy, only: Concrete_Particle_Energy, &
    particle_energy_write_legend => Concrete_Particle_Energy_write_legend, &
    particle_energy_write => Concrete_Particle_Energy_write
use procedures_visits, only: visit
use procedures_meta_factory, only: load, create, destroy

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes(2)
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Change_Counter) :: move_counters(2)
    type(Concrete_Particle_Energy) :: particles_energies(2), particles_energies_final(2)
    real(DP) :: inter_energy, inter_energy_final

    type(json_file) :: input_data
    character(len=:), allocatable :: data_field
    logical :: data_found
    integer :: energy_units(2), inter_energy_unit
    integer :: num_steps, i_step, num_moves, i_move
    integer :: move_units(2)

    call load(input_data)
    call create(environment, input_data)
    call create(mixture, input_data, environment%periodic_box)
    call create(changes, input_data, environment%periodic_box, mixture%components)
    call create(short_potentials, input_data, environment%periodic_box, mixture)

    data_field = "Monte Carlo.number of steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)
    call check_positive("test_canonical", "num_steps", num_steps)
    deallocate(data_field)

    call input_data%destroy()

    call visit(particles_energies, inter_energy, short_potentials%intras,  &
        short_potentials%inter_micro%pair, mixture%components)

    i_step = 0
    open(newunit=energy_units(1), recl=4096, file="component_1_energy.out", action="write")
    call particle_energy_write_legend(energy_units(1))
    call particle_energy_write(energy_units(1), i_step, particles_energies(1))
    open(newunit=energy_units(2), recl=4096, file="component_2_energy.out", action="write")
    call particle_energy_write_legend(energy_units(2))
    call particle_energy_write(energy_units(2), i_step, particles_energies(2))
    open(newunit=inter_energy_unit, recl=4096, file="inter_12_energy.out", action="write")
    write(inter_energy_unit, *) "# i_step    inter"
    write(inter_energy_unit, *) i_step, inter_energy

    call metropolis_factory_create(one_particle_move, environment, changes(1)%moved_positions, &
        changes(2)%moved_positions)
    call metropolis_factory_set(one_particle_move, mixture%components, short_potentials%intras, &
        short_potentials%inters)
    call metropolis_factory_set(one_particle_move, move_counters, particles_energies, inter_energy)

    open(newunit=move_units(1), recl=4096, file="component_1_move.out", action="write")
    write(move_units(1), *) "# num_hits    num_success"
    open(newunit=move_units(2), recl=4096, file="component_2_move.out", action="write")
    write(move_units(2), *) "# num_hits    num_success"

    call short_potentials%intras(1)%particles%set(short_potentials%intras(1)%pair)
    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    do i_step = 1, num_steps
        move_counters(1)%num_hits = 0; move_counters(1)%num_success = 0
        move_counters(2)%num_hits = 0; move_counters(2)%num_success = 0
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call particle_energy_write(energy_units(1), i_step, particles_energies(1))
        call particle_energy_write(energy_units(2), i_step, particles_energies(2))
        write(inter_energy_unit, *) i_step, inter_energy
        write(move_units(1), *) i_step, move_counters(1)%num_hits, move_counters(1)%num_success
        write(move_units(2), *) i_step, move_counters(2)%num_hits, move_counters(2)%num_success
    end do

    call visit(particles_energies_final, inter_energy_final, short_potentials%intras,  &
        short_potentials%inter_micro%pair, mixture%components)
    call particle_energy_write(energy_units(1), i_step-1, particles_energies_final(1))
    call particle_energy_write(energy_units(2), i_step-1, particles_energies_final(2))
    write(inter_energy_unit, *) i_step-1, inter_energy_final

    close(move_units(1))
    close(move_units(2))
    call metropolis_factory_destroy(one_particle_move)
    close(inter_energy_unit)
    close(energy_units(2))
    close(energy_units(1))
    call destroy(short_potentials)
    call destroy(changes)
    call destroy(mixture)
    call destroy(environment)

end program test_canonical
