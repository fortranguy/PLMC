module procedures_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64, iostat_end
use data_constants, only: num_dimensions
use module_data, only: test_file_exists
use procedures_errors, only: error_exit

implicit none

private
public :: increase_coordinates_size, read_coordinates

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

    subroutine read_coordinates(coordinates, filename)
        real(DP), allocatable, intent(out) :: coordinates(:, :)
        character(len=*), intent(in) :: filename

        real(DP) :: coordinate(num_dimensions)
        integer :: num_particles, i_particle
        integer :: file_unit, read_stat

        call test_file_exists(filename)
        open(newunit=file_unit, recl=4096, file=filename, status="old", action="read")
        num_particles = 0
        do
            read(file_unit, fmt=*, iostat=read_stat) coordinate
            if (read_stat == iostat_end) exit
            num_particles = num_particles + 1
        end do
        if (num_particles > 0) then
            allocate(coordinates(num_dimensions, num_particles))
            rewind(file_unit)
            do i_particle = 1, num_particles
                read(file_unit, *) coordinates(:, i_particle)
            end do
        end if
        close(file_unit)
    end subroutine read_coordinates

end module procedures_coordinates
