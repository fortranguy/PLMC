module procedures_complete_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_strings, only: max_line_length, max_word_length
use procedures_checks, only: check_file_exists, check_in_range
use types_raw_coordinates, only: Concrete_Raw_Coordinates
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector

implicit none

private
public :: complete_coordinates_read, complete_coordinates_deallocate

contains

    subroutine complete_coordinates_read(box_size, raw_coordinates, num_components, i_component, &
        filename, selector)
        real(DP), intent(inout) :: box_size(:)
        type(Concrete_Raw_Coordinates), intent(out) :: raw_coordinates
        integer, intent(in) :: num_components, i_component
        character(len=*), intent(in) :: filename
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        integer, dimension(num_components) :: nums_particles, nums_lbounds, nums_ubounds
        integer :: global_i_particle, i_particle
        integer :: coordinates_unit
        character(len=1) :: comment_character
        character(len=max_word_length) :: field

        call check_file_exists(filename)
        open(newunit=coordinates_unit, recl=max_line_length, file=filename, status="old", &
            action="read")
        read(coordinates_unit, *) comment_character, field, box_size
        read(coordinates_unit, *) comment_character, field, nums_particles
        read(coordinates_unit, *) comment_character
        call set_bounds(nums_lbounds, nums_ubounds, nums_particles)
        call allocate_coordinates(raw_coordinates, nums_particles(i_component), selector)
        i_particle = 0
        do global_i_particle = 1, sum(nums_particles)
            if (nums_lbounds(i_component) <= global_i_particle .and. global_i_particle <= &
                nums_ubounds(i_component)) then
                i_particle = i_particle + 1
                call read_coordinates(raw_coordinates, coordinates_unit, selector, i_particle)
            else
                read(coordinates_unit, *)
            end if
        end do
        close(coordinates_unit)
    end subroutine complete_coordinates_read

    pure subroutine set_bounds(nums_lbounds, nums_ubounds, nums_particles)
        integer, intent(inout) :: nums_lbounds(:), nums_ubounds(:)
        integer, intent(in) :: nums_particles(:)

        integer :: i_component

        nums_lbounds(1) = 1
        nums_ubounds(1) = nums_particles(1)
        do i_component = 2, size(nums_particles)
            nums_lbounds(i_component) = nums_ubounds(i_component - 1) + 1
            nums_ubounds(i_component) = nums_lbounds(i_component) + nums_particles(i_component) - 1
        end do
    end subroutine set_bounds

    pure subroutine allocate_coordinates(raw_coordinates, num_particles, selector)
        type(Concrete_Raw_Coordinates), intent(out) :: raw_coordinates
        integer, intent(in) :: num_particles
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        if (selector%read_positions .and. selector%read_orientations) then
            allocate(raw_coordinates%positions(num_dimensions, num_particles))
            allocate(raw_coordinates%orientations(num_dimensions, num_particles))
        else if (selector%read_positions .and. .not.selector%read_orientations) then
            allocate(raw_coordinates%positions(num_dimensions, num_particles))
            allocate(raw_coordinates%orientations(num_dimensions, 0))
        else if (.not.selector%read_positions .and. selector%read_orientations) then
            allocate(raw_coordinates%positions(num_dimensions, 0))
            allocate(raw_coordinates%orientations(num_dimensions, num_particles))
        else
            allocate(raw_coordinates%positions(num_dimensions, 0))
            allocate(raw_coordinates%orientations(num_dimensions, 0))
        end if
    end subroutine allocate_coordinates

    subroutine read_coordinates(raw_coordinates, coordinates_unit, selector, i_particle)
        type(Concrete_Raw_Coordinates), intent(inout) :: raw_coordinates
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector
        integer, intent(in) :: coordinates_unit
        integer, intent(in) :: i_particle

        integer :: i_component_dummy

        if (selector%read_positions .and. selector%read_orientations) then
            read(coordinates_unit, *) i_component_dummy, raw_coordinates%positions(:, i_particle), &
                raw_coordinates%orientations(:, i_particle)
        else if (selector%read_positions .and. .not.selector%read_orientations) then
            read(coordinates_unit, *) i_component_dummy, raw_coordinates%positions(:, i_particle)
        else if (.not.selector%read_positions .and. selector%read_orientations) then
            read(coordinates_unit, *) i_component_dummy, raw_coordinates%orientations(:, i_particle)
        end if
    end subroutine read_coordinates

    subroutine complete_coordinates_deallocate(raw_coordinates)
        type(Concrete_Raw_Coordinates), intent(inout) :: raw_coordinates

        if (allocated(raw_coordinates%orientations)) deallocate(raw_coordinates%orientations)
        if (allocated(raw_coordinates%positions)) deallocate(raw_coordinates%positions)
    end subroutine complete_coordinates_deallocate

end module procedures_complete_coordinates_reader
