module data_constants

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: num_dimensions, num_components, real_zero

    integer, parameter :: num_dimensions = 3
    integer, parameter :: num_components = 2
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: real_zero = real(2**4, DP) * epsilon(1._DP)

end module data_constants
