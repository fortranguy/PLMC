module module_particles_energy

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: Concrete_Particles_Energy_sum, operator(+), operator(-)

    type, public :: Concrete_Particles_Energy
        real(DP) :: intra = 0._DP
        real(DP) :: walls = 0._DP
    end type Concrete_Particles_Energy

    interface assignment(=)
        module procedure :: Concrete_Particles_Energy_assignment
    end interface

    interface operator(+)
        module procedure :: Concrete_Particles_Energy_addition
    end interface

    interface operator(-)
        module procedure :: Concrete_Particles_Energy_difference
    end interface

contains

    pure function Concrete_Particles_Energy_sum(particle_energy) result(particle_energy_sum)
        type(Concrete_Particles_Energy), intent(in) :: particle_energy
        real(DP) :: particle_energy_sum

        particle_energy_sum = particle_energy%intra + particle_energy%walls
    end function Concrete_Particles_Energy_sum

    pure subroutine Concrete_Particles_Energy_assignment(particle_energy_target, &
        particle_energy_value)
        type(Concrete_Particles_Energy), intent(out) :: particle_energy_target
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_value

        particle_energy_target%intra = particle_energy_value%intra
        particle_energy_target%walls = particle_energy_value%walls
    end subroutine Concrete_Particles_Energy_assignment

    function Concrete_Particles_Energy_addition(particle_energy_1, particle_energy_2) &
        result(particle_energy_addition)
        type(Concrete_Particles_Energy) :: particle_energy_addition
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_1, particle_energy_2

        particle_energy_addition%intra = particle_energy_1%intra + particle_energy_2%intra
        particle_energy_addition%walls = particle_energy_1%walls + particle_energy_2%walls
    end function Concrete_Particles_Energy_addition

    function Concrete_Particles_Energy_difference(particle_energy_1, particle_energy_2) &
        result(particle_energy_difference)
        type(Concrete_Particles_Energy) :: particle_energy_difference
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_1, particle_energy_2

        particle_energy_difference%intra = particle_energy_1%intra - particle_energy_2%intra
        particle_energy_difference%walls = particle_energy_1%walls - particle_energy_2%walls
    end function Concrete_Particles_Energy_difference

end module module_particles_energy
