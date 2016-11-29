module types_particle_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, public :: Concrete_Particle
        integer :: i = 0
        real(DP) :: position(num_dimensions) = 0._DP
        real(DP) :: orientation(num_dimensions) = 0._DP
        real(DP) :: dipole_moment(num_dimensions) = 0._DP
    end type Concrete_Particle

    type, public :: Concrete_Double_Particle
        type(Concrete_Particle) :: particles(2)
    end type Concrete_Double_Particle

end module types_particle_wrapper
