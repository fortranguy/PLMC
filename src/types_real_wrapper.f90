module types_real_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Real_Line
        real(DP), allocatable :: line(:)
    end type Real_Line

end module types_real_wrapper
