module data_constants

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none
    
    real(DP), parameter :: PI = acos(-1._DP)
    real(DP), parameter :: sigma3d = 1._DP/sqrt(3._DP)
        
end module data_constants
