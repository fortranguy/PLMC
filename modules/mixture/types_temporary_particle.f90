module types_temporary_particle

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, public :: Concrete_Temporary_Particle
        integer :: i = 0
        real(DP) :: position(num_dimensions) = 0._DP
        real(DP) :: orientation(num_dimensions) = 0._DP
        real(DP) :: dipole_moment(num_dimensions) = 0._DP
    end type Concrete_Temporary_Particle

end module types_temporary_particle
