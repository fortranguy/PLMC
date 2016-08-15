module procedures_metropolis_test

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: metropolis_test

contains

    logical function metropolis_test(acceptation_probability) result(success)
        real(DP), intent(in) :: acceptation_probability

        real(DP) :: rand

        if (acceptation_probability < 1._DP) then
            call random_number(rand)
            if (rand < acceptation_probability) then
                success = .true.
            else
                success = .false.
            end if
        else
            success = .true.
        end if
    end function metropolis_test

end module procedures_metropolis_test
