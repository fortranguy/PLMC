module procedures_random_number

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_array_size

implicit none

private
public :: random_integer, normal_random_number, markov_orientation

interface normal_random_number
    module procedure :: normal_random_number_scalar
    module procedure :: normal_random_number_vector
end interface normal_random_number

    real(DP), parameter :: sigma_3d = 1._DP / sqrt(3._DP)
    logical :: result_saved = .false.
    real(DP) :: saved_result = 0._DP

contains

    integer function random_integer(maximum)
        integer, intent(in) :: maximum

        real(DP) :: rand

        call check_positive("procedures_random_number: random_integer", "maximum", maximum)
        if (maximum == 1) then
            random_integer = maximum
        else
            call random_number(rand)
            random_integer = int(rand * real(maximum, DP)) + 1
        end if
    end function random_integer

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

    !> From SMAC, Algorithm 1.24, p.44
    subroutine markov_orientation(orientation, orientation_delta)
        real(DP), intent(inout) :: orientation(:)
        real(DP), intent(in) :: orientation_delta

        real(DP) :: rotation(num_dimensions)
        real(DP) :: amplitude, rand

        call check_array_size("procedures_random_number: markov_orientation", "orientation", &
            orientation, num_dimensions)

        call normal_random_number(rotation)
        rotation = rotation - dot_product(rotation, orientation) * orientation
        rotation = rotation / norm2(rotation)
        call random_number(rand)
        amplitude = orientation_delta * (rand - 0.5_DP)
        orientation = orientation + amplitude * rotation
        orientation = orientation / norm2(orientation)
    end subroutine markov_orientation

end module procedures_random_number
