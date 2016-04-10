module procedures_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, max_line_length, max_word_length
use procedures_checks, only: check_file_exists

implicit none

private
public :: create_positions_from_file, create_coordinates_from_file

contains

    subroutine create_positions_from_file(positions, filename)
        real(DP), dimension(:, :), allocatable, intent(out) :: positions
        character(len=*), intent(in) :: filename

        real(DP), dimension(:, :), allocatable :: orientations_dummy
        call create_coordinates_from_file(positions, orientations_dummy, filename, &
            read_orientations=.false.)
    end subroutine create_positions_from_file

    subroutine create_coordinates_from_file(positions, orientations, filename, read_orientations)
        real(DP), dimension(:, :), allocatable, intent(out) :: positions, orientations
        character(len=*), intent(in) :: filename
        logical, intent(in) :: read_orientations

        character(len=max_word_length) :: comment_caracter, number_field
        integer :: num_particles, i_particle
        integer :: file_unit

        call check_file_exists(filename)
        open(newunit=file_unit, recl=max_line_length, file=filename, status="old", action="read")
        read(file_unit, *) comment_caracter, number_field, num_particles
        read(file_unit, *) comment_caracter ! header
        allocate(positions(num_dimensions, num_particles))
        if (read_orientations) then
            allocate(orientations(num_dimensions, num_particles))
            do i_particle = 1, num_particles
                read(file_unit, *) positions(:, i_particle), orientations(:, i_particle)
            end do
        else
            do i_particle = 1, num_particles
                read(file_unit, *) positions(:, i_particle)
            end do
        end if
        close(file_unit)
    end subroutine create_coordinates_from_file

end module procedures_coordinates_reader
