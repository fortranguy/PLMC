module types_particle

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

    type, public :: Concrete_Particle
        logical :: same_type = .false.
        integer :: i = 0
        real(DP) :: diameter = 0._DP
        real(DP) :: min_diameter = 0._DP
        real(DP) :: moment_norm = 0._DP
        real(DP) :: position(num_dimensions) = 0._DP
        real(DP) :: orientation(num_dimensions) = 0._DP
    end type Concrete_Particle

end module types_particle
