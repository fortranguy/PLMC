module types_reals_line

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Reals_Line
        real(DP), allocatable :: line(:)
    end type Reals_Line

end module types_reals_line
