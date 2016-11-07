module types_real_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Real_Line
        real(DP), allocatable :: line(:)
    end type Real_Line

    type, public :: Real_Triangle
        type(Real_Line), allocatable :: triangle(:)
    end type Real_Triangle

    type, public :: Real_Triangle_Line
        type(Real_Triangle), allocatable :: line(:)
    end type Real_Triangle_Line

end module types_real_wrapper
