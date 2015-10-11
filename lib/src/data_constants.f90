module data_constants

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: num_components, num_dimensions, num_algorithms, max_word_length, max_line_length, &
    PI, real_zero

    integer, parameter :: num_components = 2
    integer, parameter :: num_dimensions = 3
    integer, parameter :: num_algorithms = 2
    integer, parameter :: max_word_length = 1024
    integer, parameter :: max_line_length = 4096
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: real_zero = real(2**4, DP) * epsilon(1._DP)

end module data_constants
