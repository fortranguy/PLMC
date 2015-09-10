module data_precisions

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public real_zero

    real(DP), parameter :: real_zero = real(2**4, DP) * epsilon(1._DP)
    real(DP), parameter :: input_tiny = real(2**2, DP) * epsilon(1._DP) ! to move?
    real(DP), parameter :: consistency_tiny = real(2**13, DP) * epsilon(1._DP) ! to move?
    real(DP), parameter :: field_consistency_tiny = real(2**17, DP) * epsilon(1._DP) ! to move?

end module data_precisions
