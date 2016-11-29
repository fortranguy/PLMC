module data_constants

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: num_dimensions, PI, real_zero

    integer, parameter :: num_dimensions = 3
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: real_zero = real(2**5, DP) * epsilon(1._DP)

end module data_constants
