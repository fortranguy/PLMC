module procedures_normal_random_number

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: normal_random_number

interface normal_random_number
    module procedure :: normal_random_number_scalar
    module procedure :: normal_random_number_vector
end interface normal_random_number

    real(DP), parameter :: sigma_3d = 1._DP / sqrt(3._DP)
    logical :: result_saved = .false.
    real(DP) :: saved_result

contains

    !> From SMAC, Algorithm 1.19, p.39
    subroutine normal_random_number_scalar(gauss)
        real(DP), intent(out) :: gauss

        real(DP), dimension(2) :: xy, rand_2d
        real(DP) :: norm_squared, factor

        if (.not.result_saved) then
            do
                call random_number(rand_2d)
                xy = 2._DP*rand_2d - 1._DP
                norm_squared = sum(xy**2)
                if (0._DP < norm_squared .and. norm_squared <= 1._DP) exit
            end do
            factor = sqrt(-2._DP * log(norm_squared) / norm_squared)
            gauss = sigma_3d * factor * xy(1)
            saved_result = sigma_3d * factor * xy(2)
            result_saved = .true.
        else
            gauss = saved_result
            result_saved = .false.
        end if
    end subroutine normal_random_number_scalar

    subroutine normal_random_number_vector(gauss)
        real(DP), intent(inout) :: gauss(:)

        integer :: i_dimension

        do i_dimension = 1, size(gauss)
            call normal_random_number(gauss(i_dimension))
        end do
    end subroutine normal_random_number_vector

end module procedures_normal_random_number
