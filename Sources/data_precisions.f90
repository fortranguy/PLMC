module data_precisions

use, intrinsic :: iso_fortran_env, only: REAL64

implicit none

    integer, parameter :: DP = REAL64 !< double precision
    real(DP), parameter :: real_zero = real(2**4, DP) * epsilon(1._DP)
    real(DP), parameter :: io_tiny = real(2**2, DP) * epsilon(1._DP)
    real(DP), parameter :: consist_tiny = real(2**13, DP) * epsilon(1._DP)

end module data_precisions
