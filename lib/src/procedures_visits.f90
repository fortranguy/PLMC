module procedures_visits

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use class_walls_potential, only: Abstract_Walls_Potential
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use types_particles_wrapper, only: Particles_Wrapper
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use types_short_potential_wrapper, only: Short_Potential_Wrapper
use module_particles_energy, only: Concrete_Particles_Energy

implicit none

private
public :: visit

interface visit
    module procedure :: visit_mixture
    module procedure :: visit_particles
    module procedure :: visit_particles_walls
end interface visit

contains

    subroutine visit_mixture_walls(particles_energies, walls_potential, components, wall_pairs)
        type(Concrete_Particles_Energy), intent(inout) :: particles_energies(2)
        class(Abstract_Walls_Potential), intent(inout) :: walls_potential
        type(Particles_Wrapper), intent(in) :: components(2)
        class(Abstract_Pair_Potential), intent(in) :: wall_pairs(2)

        logical :: overlap

        call walls_potential%set(wall_pairs(1))
        call visit(overlap, particles_energies(1)%walls, walls_potential, components(1)%positions)
        if (overlap) call error_exit("walls - components(1) overlap")
        call walls_potential%set(wall_pairs(2))
        call visit(overlap, particles_energies(2)%walls, walls_potential, components(2)%positions)
        if (overlap) call error_exit("walls - components(1) overlap")
    end subroutine visit_mixture_walls

    pure subroutine visit_particles_walls(overlap, energy, walls_potential, particles_positions)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        class(Abstract_Particles_Positions), intent(in) :: particles_positions

        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy_i = 0._DP
        do i_particle = 1, particles_positions%get_num()
            call walls_potential%visit(overlap, energy_i, particles_positions%get(i_particle))
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine visit_particles_walls

    subroutine visit_mixture(particles_energies, inter_energy, intras, inter_pair, &
        components)
        type(Concrete_Particles_Energy), intent(inout) :: particles_energies(2)
        real(DP), intent(out) :: inter_energy
        type(Short_Potential_Wrapper), intent(inout) :: intras(2)
        class(Abstract_Pair_Potential), intent(in) :: inter_pair
        type(Particles_Wrapper), intent(in) :: components(2)

        logical :: overlap
        real(DP) :: inter_energy_1, inter_energy_2

        call intras(1)%particles%set(intras(1)%pair)
        call visit(overlap, particles_energies(1)%intra, intras(1)%particles, &
            components(1)%positions, same_type=.true.)
        if (overlap) call error_exit("intra 1 overlap")
        call intras(2)%particles%set(intras(2)%pair)
        call visit(overlap, particles_energies(2)%intra, intras(2)%particles, &
            components(2)%positions, same_type=.true.)
        if (overlap) call error_exit("intra 2 overlap")
        call intras(1)%particles%set(inter_pair)
        call visit(overlap, inter_energy_1, intras(1)%particles, components(2)%positions, &
            same_type=.false.)
        if (overlap) call error_exit("inter intras(1) overlap")
        call intras(2)%particles%set(inter_pair)
        call visit(overlap, inter_energy_2, intras(2)%particles, components(1)%positions, &
        same_type=.false.)
        if (overlap) call error_exit("inter intras(2) overlap")
        inter_energy = inter_energy_1 + inter_energy_2
    end subroutine visit_mixture

    pure subroutine visit_particles(overlap, energy, particles_potential, &
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
            if (overlap) return
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit_particles

end module procedures_visits
