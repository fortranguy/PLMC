module procedures_random

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: random_integer

contains

    function random_integer(maximum)
        integer, intent(in) :: maximum
        integer :: random_integer

        real(DP) :: rand
        call random_number(rand)
        random_integer = int(rand * real(maximum, DP)) + 1
    end function random_integer

end module procedures_random
