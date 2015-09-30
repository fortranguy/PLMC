module procedures_visits

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_potential, only: Abstract_Particles_Potential

implicit none

private
public :: visit

interface visit
    module procedure :: particles_potential_visit
end interface visit

contains

    pure subroutine particles_potential_visit(overlap, energy, particles_potential, &
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
    end subroutine particles_potential_visit

end module procedures_visits
