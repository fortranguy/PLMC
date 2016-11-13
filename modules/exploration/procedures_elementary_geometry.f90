module procedures_elementary_geometry

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI

implicit none

private
public :: sphere_surface

contains

    !> \[ S(r) = 4\pi r^2 \]
    pure real(DP) function sphere_surface(radius)
        real(DP), intent(in) :: radius

        sphere_surface = 4._DP*PI * radius**2
    end function sphere_surface

end module procedures_elementary_geometry
