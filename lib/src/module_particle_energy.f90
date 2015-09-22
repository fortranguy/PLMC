module module_particle_energy

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: Concrete_Particle_Energy_sum, operator(-)

    type, public :: Concrete_Particle_Energy
        real(DP) :: intra = 0._DP
        real(DP) :: inter = 0._DP
    end type Concrete_Particle_Energy

    interface assignment(=)
        module procedure Concrete_Particle_Energy_assignment
    end interface

    interface operator(-)
        module procedure Concrete_Particle_Energy_difference
    end interface

contains

    pure function Concrete_Particle_Energy_sum(particle_energy) result(particle_energy_sum)
        type(Concrete_Particle_Energy), intent(in) :: particle_energy
        real(DP) :: particle_energy_sum

        particle_energy_sum = particle_energy%intra + particle_energy%inter
    end function Concrete_Particle_Energy_sum

    pure subroutine Concrete_Particle_Energy_assignment(particle_energy_target, &
        particle_energy_value)
        type(Concrete_Particle_Energy), intent(out) :: particle_energy_target
        type(Concrete_Particle_Energy), intent(in) :: particle_energy_value

        particle_energy_target%intra = particle_energy_value%intra
        particle_energy_target%inter = particle_energy_value%inter
    end subroutine Concrete_Particle_Energy_assignment

    function Concrete_Particle_Energy_difference(particle_energy_1, particle_energy_2) &
        result(particle_energy_difference)
        type(Concrete_Particle_Energy) :: particle_energy_difference
        type(Concrete_Particle_Energy), intent(in) :: particle_energy_1, particle_energy_2

        particle_energy_difference%intra = particle_energy_1%intra - particle_energy_2%intra
        particle_energy_difference%inter = particle_energy_1%inter - particle_energy_2%inter
    end function Concrete_Particle_Energy_difference

end module module_particle_energy