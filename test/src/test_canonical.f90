module procedures_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_potential, only: Abstract_Particles_Potential

implicit none

private
public :: visit

contains

    subroutine visit(overlap, energy, particles_potential, particles_positions, same_type)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Particles_Potential), intent(in) :: particles_potential
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        logical, intent(in) :: same_type

        type(Concrete_Particle) :: particle
        real(DP) :: energy_i
        integer :: i_particle

        particle%same_type = same_type
        overlap = .false.
        energy = 0._DP
        do i_particle = 1, particles_positions%get_num()
            particle%i = i_particle
            particle%position = particles_positions%get(particle%i)
            call particles_potential%visit(overlap, energy_i, particle)
            if (overlap) return
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit

end module procedures_canonical

program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use types_particles_wrapper, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy
use types_short_potential, only: Mixture_Short_Potentials_Wrapper
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use class_one_particle_move, only: Abstract_One_Particle_Move
use procedures_metropolis_factory, only: metropolis_factory_create, metropolis_factory_destroy
use types_change_counter, only: Concrete_Change_Counter
use module_particle_energy, only: Concrete_Particle_Energy, &
    particle_energy_write_legend => Concrete_Particle_Energy_write_legend, &
    particle_energy_write => Concrete_Particle_Energy_write
use procedures_canonical, only: visit

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes_1, changes_2
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Change_Counter) :: move_counters(2)
    type(Concrete_Particle_Energy) :: particles_energies(2)
    real(DP) :: inter_energy, inter_energy_1, inter_energy_2, particles_energy_1_intra
    logical :: overlap

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    integer :: energy_units(2), inter_energy_unit
    integer :: num_steps, i_step, num_moves, i_move
    integer :: move_units(2)

    call json_initialize()
    data_filename = "canonical.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call environment_factory_create(environment, input_data, "Environment.")

    call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
        environment%periodic_box)
    call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
        environment%periodic_box)
    call particles_factory_create(mixture%inter_diameter, mixture%components(1)%diameter, &
        mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")

    call changes_factory_create(changes_1, input_data, "Changes.Component 1.", &
        environment%periodic_box, mixture%components(1))
    call changes_factory_create(changes_2, input_data, "Changes.Component 2.", &
        environment%periodic_box, mixture%components(2))

    call short_potential_factory_create(short_potentials%intras(1), input_data, &
        "Short Potentials.Component 1.", environment%periodic_box, mixture%components(1))
    call short_potential_factory_create(short_potentials%intras(2), input_data, &
        "Short Potentials.Component 2.", environment%periodic_box, mixture%components(2))
    call short_potential_factory_create(short_potentials%inter_micro, input_data, &
        "Short Potentials.Inter 12.", mixture%inter_diameter)
    call short_potential_factory_create(short_potentials%inters(1), short_potentials%inter_micro, &
        input_data, "Short Potentials.Inter 12.", environment%periodic_box, &
        mixture%components(1)%positions)
    call short_potential_factory_create(short_potentials%inters(2), short_potentials%inter_micro, &
        input_data, "Short Potentials.Inter 12.", environment%periodic_box, &
        mixture%components(2)%positions)

    data_field = "Monte Carlo.number of steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)
    call check_positive("test_canonical", "num_steps", num_steps)
    deallocate(data_field)

    call input_data%destroy()

    call short_potentials%intras(1)%particles%set(short_potentials%intras(1)%pair)
    call visit(overlap, particles_energies(1)%intra, short_potentials%intras(1)%particles, &
        mixture%components(1)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(1) overlap")
    call short_potentials%intras(2)%particles%set(short_potentials%intras(2)%pair)
    call visit(overlap, particles_energies(2)%intra, short_potentials%intras(2)%particles, &
        mixture%components(2)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(2) overlap")
    call short_potentials%intras(1)%particles%set(short_potentials%inter_micro%pair)
    call visit(overlap, inter_energy_1, short_potentials%intras(1)%particles, &
        mixture%components(2)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(1) overlap")
    call short_potentials%intras(2)%particles%set(short_potentials%inter_micro%pair)
    call visit(overlap, inter_energy_2, short_potentials%intras(2)%particles, &
        mixture%components(1)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(2) overlap")
    inter_energy = inter_energy_1 + inter_energy_2

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

    call metropolis_factory_create(one_particle_move, environment, changes_1%moved_positions, &
        changes_2%moved_positions)
    call one_particle_move%set_candidate(1, mixture%components(1)%positions)
    call one_particle_move%set_candidate(1, short_potentials%intras(1)%cells, &
        short_potentials%inters(1)%cells)
    call one_particle_move%set_candidate(2, mixture%components(2)%positions)
    call one_particle_move%set_candidate(2, short_potentials%intras(2)%cells, &
        short_potentials%inters(2)%cells)
    call one_particle_move%set_candidates_observables(move_counters, particles_energies, &
        inter_energy)

    open(newunit=move_units(1), recl=4096, file="component_1_move.out", action="write")
    write(move_units(1), *) "# num_hits    num_success"
    open(newunit=move_units(2), recl=4096, file="component_2_move.out", action="write")
    write(move_units(2), *) "# num_hits    num_success"

    call short_potentials%intras(1)%particles%set(short_potentials%intras(1)%pair)
    num_moves = mixture%components(1)%positions%get_num() + &
        mixture%components(2)%positions%get_num()
    do i_step = 1, num_steps
        do i_move = 1, num_moves
            call one_particle_move%try()
        end do
        call particle_energy_write(energy_units(1), i_step, particles_energies(1))
        call particle_energy_write(energy_units(2), i_step, particles_energies(2))
        write(inter_energy_unit, *) i_step, inter_energy
        write(move_units(1), *) i_step, move_counters(1)%num_hits, move_counters(1)%num_success
        write(move_units(1), *) i_step, move_counters(2)%num_hits, move_counters(2)%num_success
    end do

    call short_potentials%intras(1)%particles%set(short_potentials%intras(1)%pair)
    call visit(overlap, particles_energies(1)%intra, short_potentials%intras(1)%particles, &
        mixture%components(1)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(1) overlap")
    call short_potentials%intras(2)%particles%set(short_potentials%intras(2)%pair)
    call visit(overlap, particles_energies(2)%intra, short_potentials%intras(2)%particles, &
        mixture%components(2)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(2) overlap")
    call short_potentials%intras(1)%particles%set(short_potentials%inter_micro%pair)
    call visit(overlap, inter_energy_1, short_potentials%intras(1)%particles, &
        mixture%components(2)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(1) overlap")
    call short_potentials%intras(2)%particles%set(short_potentials%inter_micro%pair)
    call visit(overlap, inter_energy_2, short_potentials%intras(2)%particles, &
        mixture%components(1)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(2) overlap")
    inter_energy = inter_energy_1 + inter_energy_2
    call particle_energy_write(energy_units(1), i_step, particles_energies(1))
    call particle_energy_write(energy_units(2), i_step, particles_energies(2))
    write(inter_energy_unit, *) i_step, inter_energy

    close(move_units(1))
    close(move_units(2))
    call metropolis_factory_destroy(one_particle_move)
    close(inter_energy_unit)
    close(energy_units(2))
    close(energy_units(1))
    call short_potential_factory_destroy(short_potentials%inters(2))
    call short_potential_factory_destroy(short_potentials%inters(1))
    call short_potential_factory_destroy(short_potentials%inter_micro)
    call short_potential_factory_destroy(short_potentials%intras(2))
    call short_potential_factory_destroy(short_potentials%intras(1))
    call changes_factory_destroy(changes_1)
    call changes_factory_destroy(changes_2)
    call particles_factory_destroy(mixture%inter_diameter)
    call particles_factory_destroy(mixture%components(2))
    call particles_factory_destroy(mixture%components(1))
    call environment_factory_destroy(environment)

end program test_canonical
