module module_particles_energy

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: Concrete_Particles_Energy_sum, Concrete_Inter_Energy_sum, &
    operator(+), operator(-)

    type, public :: Concrete_Particles_Energy
        real(DP) :: short = 0._DP
        real(DP) :: walls = 0._DP
        real(DP) :: long = 0._DP
    end type Concrete_Particles_Energy

    type, public :: Concrete_Inter_Energy
        real(DP) :: short = 0._DP
        real(DP) :: long = 0._DP
    end type Concrete_Inter_Energy

    interface assignment(=)
        module procedure :: Concrete_Particles_Energy_assignment
        module procedure :: Concrete_Inter_Energy_assignment
    end interface

    interface operator(+)
        module procedure :: Concrete_Particles_Energy_addition
        module procedure :: Concrete_Inter_Energy_addition
    end interface

    interface operator(-)
        module procedure :: Concrete_Particles_Energy_difference
        module procedure :: Concrete_Inter_Energy_difference
    end interface

contains

!implementation Concrete_Particles_Energy

    pure function Concrete_Particles_Energy_sum(particle_energy) result(particle_energy_sum)
        type(Concrete_Particles_Energy), intent(in) :: particle_energy
        real(DP) :: particle_energy_sum

        particle_energy_sum = particle_energy%short + particle_energy%walls + particle_energy%long
    end function Concrete_Particles_Energy_sum

    pure subroutine Concrete_Particles_Energy_assignment(particle_energy_target, &
        particle_energy_value)
        type(Concrete_Particles_Energy), intent(out) :: particle_energy_target
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_value

        particle_energy_target%short = particle_energy_value%short
        particle_energy_target%walls = particle_energy_value%walls
        particle_energy_target%long = particle_energy_value%long
    end subroutine Concrete_Particles_Energy_assignment

    function Concrete_Particles_Energy_addition(particle_energy_1, particle_energy_2) &
        result(particle_energy_addition)
        type(Concrete_Particles_Energy) :: particle_energy_addition
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_1, particle_energy_2

        particle_energy_addition%short = particle_energy_1%short + particle_energy_2%short
        particle_energy_addition%walls = particle_energy_1%walls + particle_energy_2%walls
        particle_energy_addition%long = particle_energy_1%long + particle_energy_2%long
    end function Concrete_Particles_Energy_addition

    function Concrete_Particles_Energy_difference(particle_energy_1, particle_energy_2) &
        result(particle_energy_difference)
        type(Concrete_Particles_Energy) :: particle_energy_difference
        type(Concrete_Particles_Energy), intent(in) :: particle_energy_1, particle_energy_2

        particle_energy_difference%short = particle_energy_1%short - particle_energy_2%short
        particle_energy_difference%walls = particle_energy_1%walls - particle_energy_2%walls
        particle_energy_difference%long = particle_energy_1%long - particle_energy_2%long
    end function Concrete_Particles_Energy_difference

!end implementation Concrete_Particles_Energy

!implementation Concrete_Inter_Energy

    pure function Concrete_Inter_Energy_sum(inter_energy) result(inter_energy_sum)
        type(Concrete_Inter_Energy), intent(in) :: inter_energy
        real(DP) :: inter_energy_sum

        inter_energy_sum = inter_energy%short + inter_energy%long
    end function Concrete_Inter_Energy_sum

    pure subroutine Concrete_Inter_Energy_assignment(inter_energy_target, &
        inter_energy_value)
        type(Concrete_Inter_Energy), intent(out) :: inter_energy_target
        type(Concrete_Inter_Energy), intent(in) :: inter_energy_value

        inter_energy_target%short = inter_energy_value%short
        inter_energy_target%long = inter_energy_value%long
    end subroutine Concrete_Inter_Energy_assignment

    function Concrete_Inter_Energy_addition(inter_energy_1, inter_energy_2) &
        result(inter_energy_addition)
        type(Concrete_Inter_Energy) :: inter_energy_addition
        type(Concrete_Inter_Energy), intent(in) :: inter_energy_1, inter_energy_2

        inter_energy_addition%short = inter_energy_1%short + inter_energy_2%short
        inter_energy_addition%long = inter_energy_1%long + inter_energy_2%long
    end function Concrete_Inter_Energy_addition

    function Concrete_Inter_Energy_difference(inter_energy_1, inter_energy_2) &
        result(inter_energy_difference)
        type(Concrete_Inter_Energy) :: inter_energy_difference
        type(Concrete_Inter_Energy), intent(in) :: inter_energy_1, inter_energy_2

        inter_energy_difference%short = inter_energy_1%short - inter_energy_2%short
        inter_energy_difference%long = inter_energy_1%long - inter_energy_2%long
    end function Concrete_Inter_Energy_difference

!end implementation Concrete_Inter_Energy

end module module_particles_energy
