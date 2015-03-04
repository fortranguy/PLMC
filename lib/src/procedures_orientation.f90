module procedures_orientation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

public random_orientation, markov_orientation

    real(DP), parameter :: sigma3d = 1._DP / sqrt(3._DP)

contains

    !> From SMAC, Algorithm 1.23, p. 43
    function random_orientation()
        real(DP), dimension(num_dimensions) :: random_orientation
        
        integer :: i_dimension
        
        do i_dimension = 1, num_dimensions
            random_orientation(i_dimension) = gauss()
        end do
        random_orientation = random_orientation / norm2(random_orientation)
    end function random_orientation

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
    
    !> From SMAC, Algorithm 1.24, p.44
    subroutine markov_orientation(orientation, orientation_delta)
        real(DP), intent(inout) :: orientation(num_dimensions)
        real(DP), intent(in) :: orientation_delta
        
        real(DP) :: rotation(num_dimensions)
        real(DP) :: rotation_dot_orientation
        real(DP) :: amplitude, random
        integer :: i_dimension
        
        do i_dimension = 1, num_dimensions
            rotation(i_dimension) = gauss()
        end do
        
        rotation_dot_orientation = dot_product(rotation, orientation)
        rotation = rotation - rotation_dot_orientation * orientation
        rotation = rotation / norm2(rotation)
        
        call random_number(random)
        amplitude = orientation_delta * (random - 0.5_DP)
        
        orientation = orientation + amplitude * rotation
        orientation = orientation / norm2(orientation)
    end subroutine markov_orientation

end module procedures_orientation
