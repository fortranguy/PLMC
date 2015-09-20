module types_metropolis_algorithms

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Particle_Energy
        real(DP) :: intra
        real(DP) :: inter
        real(DP) :: sum
    end type Concrete_Particle_Energy

end module types_metropolis_algorithms
