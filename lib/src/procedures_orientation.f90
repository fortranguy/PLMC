module procedures_orientation

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

public gauss

    real(DP), parameter :: sigma3d = 1._DP / sqrt(3._DP)

contains

    !> From SMAC, Algorithm 1.19, p.39
    function gauss()
        real(DP) :: gauss

        real(DP) :: x, y
        real(DP) :: norm_sqr, factor
        real(DP), save :: gauss_save
        integer, save :: switch = 0

        if (switch == 0) then
            do
                call random_number(x)
                x = 2._DP*x - 1._DP
                call random_number(y)
                y = 2._DP*y - 1._DP

                norm_sqr = x*x + y*y

                if (norm_sqr <= 1._DP .and. norm_sqr > 0._DP) exit
            end do
            factor = sqrt(-2._DP * log(norm_sqr) / norm_sqr)
            gauss_save = sigma3d * factor * x
            gauss = sigma3d * factor * y
            switch = 1
        else
            gauss = gauss_save
            switch = 0
        end if
    end function gauss

end module procedures_orientation