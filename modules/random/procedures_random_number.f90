module procedures_random_number

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive, check_array_size
use procedures_normal_random_number, only: normal_random_number

implicit none

private
public :: random_integer, markov_orientation

contains

    function random_integer(maximum)
        integer, intent(in) :: maximum
        integer :: random_integer

        real(DP) :: rand

        call check_positive("procedures_random_number: random_integer", "maximum", maximum)
        if (maximum == 1) then
            random_integer = maximum
        else
            call random_number(rand)
            random_integer = int(rand * real(maximum, DP)) + 1
        end if
    end function random_integer

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
