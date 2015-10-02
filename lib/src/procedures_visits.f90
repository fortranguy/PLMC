module procedures_visits

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use class_walls_potential, only: Abstract_Walls_Potential
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use types_particles_wrapper, only: Particles_Wrapper
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use types_short_potential_wrapper, only: Short_Potential_Wrapper, Wall_Short_Potential_Wrapper
use module_particles_energy, only: Concrete_Particles_Energy

implicit none

private
public :: visit

interface visit
    module procedure :: visit_mixture
    module procedure :: visit_particles
    module procedure :: visit_particles_walls
    module procedure :: visit_mixture_walls
end interface visit

contains

    subroutine visit_mixture_walls(particles_energies, walls_potential, components, walls)
        type(Concrete_Particles_Energy), intent(inout) :: particles_energies(2)
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        type(Particles_Wrapper), intent(in) :: components(2)
        class(Wall_Short_Potential_Wrapper), intent(in) :: walls(2)

        logical :: overlap

        call visit(overlap, particles_energies(1)%walls, walls_potential, components(1)%positions, &
            walls(1)%pair)
        if (overlap) call error_exit("walls - components(1) overlap")
        call visit(overlap, particles_energies(2)%walls, walls_potential, components(2)%positions, &
            walls(2)%pair)
        if (overlap) call error_exit("walls - components(2) overlap")
    end subroutine visit_mixture_walls

    pure subroutine visit_particles_walls(overlap, energy, walls_potential, particles_positions, &
        wall_pair)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Pair_Potential), intent(in) :: wall_pair

        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, particles_positions%get_num()
            call walls_potential%visit(overlap, energy_i, particles_positions%get(i_particle), &
                wall_pair)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine visit_particles_walls

    subroutine visit_mixture(particles_energies, inter_energy, intras, inter_pair, &
        components)
        type(Concrete_Particles_Energy), intent(inout) :: particles_energies(2)
        real(DP), intent(out) :: inter_energy
        type(Short_Potential_Wrapper), intent(in) :: intras(2)
        class(Abstract_Pair_Potential), intent(in) :: inter_pair
        type(Particles_Wrapper), intent(in) :: components(2)

        logical :: overlap
        real(DP) :: inter_energy_1, inter_energy_2

        call visit(overlap, particles_energies(1)%intra, intras(1)%particles, &
            components(1)%positions, intras(1)%pair, same_type=.true.)
        if (overlap) call error_exit("intra 1 overlap")
        call visit(overlap, particles_energies(2)%intra, intras(2)%particles, &
            components(2)%positions, intras(2)%pair, same_type=.true.)
        if (overlap) call error_exit("intra 2 overlap")
        call visit(overlap, inter_energy_1, intras(1)%particles, components(2)%positions, &
            inter_pair, same_type=.false.)
        if (overlap) call error_exit("inter intras(1) overlap")
        call visit(overlap, inter_energy_2, intras(2)%particles, components(1)%positions, &
            inter_pair, same_type=.false.)
        if (overlap) call error_exit("inter intras(2) overlap")
        inter_energy = inter_energy_1 + inter_energy_2
    end subroutine visit_mixture

    pure subroutine visit_particles(overlap, energy, particles_potential, &
        particles_positions, pair_potential, same_type)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Particles_Potential), intent(in) :: particles_potential
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
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
            call particles_potential%visit(overlap, energy_i, particle, pair_potential)
            if (overlap) return
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit_particles

end module procedures_visits
