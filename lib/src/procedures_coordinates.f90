module procedures_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions

implicit none

private

public increase_coordinates_size

    integer, parameter :: increase_factor = 2

contains

    pure subroutine increase_coordinates_size(coordinates)
        real(DP), allocatable, intent(inout) :: coordinates(:, :)
        
        real(DP), allocatable :: coordinates_tmp(:, :)
        integer :: size1, size2
        
        size1 = size(coordinates, 1)
        size2 = size(coordinates, 2)
        allocate(coordinates_tmp(size1, size2))
        coordinates_tmp = coordinates
        deallocate(coordinates)
        
        allocate(coordinates(size1, increase_factor * size2))
        coordinates(:, 1:size2) = coordinates_tmp
        deallocate(coordinates_tmp)
    end subroutine increase_coordinates_size

end module procedures_coordinates
