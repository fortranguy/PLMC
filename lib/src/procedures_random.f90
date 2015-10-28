module procedures_random

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use class_normal_random_number, only: Normal_Random_Number

implicit none

private
public :: random_integer, random_orientation, markov_orientation

contains

    function random_integer(maximum)
        integer, intent(in) :: maximum
        integer :: random_integer

        real(DP) :: rand

        call check_positive("procedures_random: random_integer", "maximum", maximum)
        if (maximum == 1) then
            random_integer = maximum
        else
            call random_number(rand)
            random_integer = int(rand * real(maximum, DP)) + 1
        end if
    end function random_integer

    !> From SMAC, Algorithm 1.23, p. 43
    function random_orientation()
        real(DP), dimension(num_dimensions) :: random_orientation

        integer :: i_dimension
        type(Normal_Random_Number) :: gauss

        do i_dimension = 1, num_dimensions
            random_orientation(i_dimension) = gauss%get()
        end do
        random_orientation = random_orientation / norm2(random_orientation)
    end function random_orientation

    !> From SMAC, Algorithm 1.24, p.44
    subroutine markov_orientation(orientation, orientation_delta)
        real(DP), intent(inout) :: orientation(:)
        real(DP), intent(in) :: orientation_delta

        real(DP) :: rotation(num_dimensions)
        type(Normal_Random_Number) :: gauss
        real(DP) :: amplitude, rand
        integer :: i_dimension

        do i_dimension = 1, num_dimensions
            rotation(i_dimension) = gauss%get()
        end do

        rotation = rotation - dot_product(rotation, orientation) * orientation
        rotation = rotation / norm2(rotation)

        call random_number(rand)
        amplitude = orientation_delta * (rand - 0.5_DP)

        orientation = orientation + amplitude * rotation
        orientation = orientation / norm2(orientation)
    end subroutine markov_orientation

end module procedures_random
