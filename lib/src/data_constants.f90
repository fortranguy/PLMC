module data_constants

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public PI
    
    real(DP), parameter :: PI = acos(-1._DP)
        
end module data_constants
