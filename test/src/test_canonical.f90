module procedures_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_potential, only: Abstract_Particles_Potential

implicit none

private
public :: calculate_particles_energy

contains

    subroutine calculate_particles_energy(overlap, energy, particles_potential, &
        particles_positions, same_type)
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
            if (overlap) exit
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine calculate_particles_energy

end module procedures_canonical

program test_canonical

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use types_box, only: Box_Wrapper
use procedures_box_factory, only: box_factory_create, box_factory_destroy
use types_particles, only: Mixture_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy
use types_short_potential, only: Mixture_Short_Potentials_Wrapper
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use class_one_particle_move, only: Abstract_One_Particle_Move
use procedures_metropolis_factory, only: metropolis_factory_create, metropolis_factory_destroy
use types_change_counter, only: Concrete_Change_Counter
use module_particle_energy, only: Concrete_Particle_Energy
use procedures_canonical, only: calculate_particles_energy

implicit none

    type(Box_Wrapper) :: box
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes_1, changes_2
    type(Mixture_Short_Potentials_Wrapper) :: short_potentials
    class(Abstract_One_Particle_Move), allocatable :: one_particle_move
    type(Concrete_Change_Counter) :: move_counters(2)
    type(Concrete_Particle_Energy) :: particles_energies(2)
    real(DP) :: inter_energy, inter_energy_1, inter_energy_2

    logical :: overlap

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "canonical.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory_create(box, input_data, "Box.")

    call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
        box%periodic_box)
    call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
        box%periodic_box)
    call particles_factory_create(mixture%inter_diameter, mixture%components(1)%diameter, &
        mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")

    call changes_factory_create(changes_1, input_data, "Changes.Component 1.", &
        mixture%components(1))
    call changes_factory_create(changes_2, input_data, "Changes.Component 2.", &
        mixture%components(2))

    call short_potential_factory_create(short_potentials%intras(1), input_data, &
        "Short Potentials.Component 1.", box%periodic_box, mixture%components(1))
    call short_potentials%intras(1)%particles%set(short_potentials%intras(1)%pair)
    call calculate_particles_energy(overlap, particles_energies(1)%intra, &
        short_potentials%intras(1)%particles, mixture%components(1)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(1) overlap")

    call short_potential_factory_create(short_potentials%intras(2), input_data, &
        "Short Potentials.Component 2.", box%periodic_box, mixture%components(2))
    call short_potentials%intras(2)%particles%set(short_potentials%intras(2)%pair)
    call calculate_particles_energy(overlap, particles_energies(2)%intra, &
        short_potentials%intras(2)%particles, mixture%components(2)%positions, same_type=.true.)
    if (overlap) call error_exit("short_potentials%intras(2) overlap")

    call short_potential_factory_create(short_potentials%inter_micro, input_data, &
        "Short Potentials.Inter 12.", mixture%inter_diameter)
    call short_potential_factory_create(short_potentials%inters(1), short_potentials%inter_micro, &
        input_data, "Short Potentials.Inter 12.", box%periodic_box, mixture%components(1)%positions)
    call short_potentials%intras(1)%particles%set(short_potentials%inter_micro%pair)
    call calculate_particles_energy(overlap, inter_energy_1, short_potentials%intras(1)%particles, &
        mixture%components(2)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(1) overlap")
    call short_potential_factory_create(short_potentials%inters(2), short_potentials%inter_micro, &
        input_data, "Short Potentials.Inter 12.", box%periodic_box, mixture%components(2)%positions)
    call short_potentials%intras(2)%particles%set(short_potentials%inter_micro%pair)
    call calculate_particles_energy(overlap, inter_energy_2, short_potentials%intras(2)%particles, &
        mixture%components(1)%positions, same_type=.false.)
    if (overlap) call error_exit("inter short_potentials%intras(2) overlap")
    inter_energy = inter_energy_1 + inter_energy_2

    call input_data%destroy()

    call metropolis_factory_create(one_particle_move, box, changes_1%moved_positions, &
        changes_2%moved_positions)
    call one_particle_move%set_candidate(1, mixture%components(1)%positions)
    call one_particle_move%set_candidate(1, short_potentials%intras(1)%cells, &
        short_potentials%inters(1)%cells)
    call one_particle_move%set_candidate(2, mixture%components(2)%positions)
    call one_particle_move%set_candidate(2, short_potentials%intras(2)%cells, &
        short_potentials%inters(2)%cells)
    inter_energy = 0._DP
    call one_particle_move%set_candidates_observables(move_counters, particles_energies, &
        inter_energy)

    call one_particle_move%try()

    call metropolis_factory_destroy(one_particle_move)
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
    call box_factory_destroy(box)

end program test_canonical
