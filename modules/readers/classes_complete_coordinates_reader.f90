module classes_complete_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_strings, only: max_line_length, max_word_length
use procedures_errors, only: error_exit
use procedures_checks, only: check_file_exists
use types_string_wrapper, only: String_Wrapper
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use classes_component_coordinates_reader, only: Component_Coordinates_Reader_wrapper
use procedures_component_coordinates_reader_factory, only: component_coordinates_reader_destroy => &
    destroy

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Reader
    private
        class(Abstract_Periodic_Box), pointer :: periodic_boxes(:) => null()
        class(Abstract_Box_Size_Checker), pointer :: boxes_size_checker(:) => null()
        type(Component_Coordinates_Reader_wrapper), allocatable :: components_coordinates(:, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: read => Abstract_read
    end type Abstract_Complete_Coordinates_Reader

    type, extends(Abstract_Complete_Coordinates_Reader), public :: &
        Concrete_Complete_Coordinates_Reader
    end type Concrete_Complete_Coordinates_Reader

contains

    subroutine Abstract_construct(this, periodic_boxes, boxes_size_checker, components_coordinates)
        class(Abstract_Complete_Coordinates_Reader), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_boxes(:)
        class(Abstract_Box_Size_Checker), target, intent(in) :: boxes_size_checker(:)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_coordinates(:, :)

        this%periodic_boxes => periodic_boxes
        this%boxes_size_checker => boxes_size_checker
        allocate(this%components_coordinates, source=components_coordinates)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Reader), intent(out) :: this

        call component_coordinates_reader_destroy(this%components_coordinates)
        this%boxes_size_checker => null()
        this%periodic_boxes => null()
    end subroutine Abstract_destroy

    !> @warning gfortran requires intent(inout) even though this%periodic_boxes is
    !> implicitly mutable with ifort.
    subroutine Abstract_read(this, coordinates)
        class(Abstract_Complete_Coordinates_Reader), intent(inout) :: this
        type(String_Wrapper), intent(in) :: coordinates(:)

        real(DP) :: box_size(num_dimensions)
        integer :: i_box, i_component
        integer :: nums_particles(size(this%components_coordinates, 1))
        character(len=1) :: comment_character
        character(len=max_word_length) :: field
        integer :: coordinates_unit

        if (size(coordinates) /= size(this%components_coordinates)) then
            call error_exit("Abstract_Complete_Coordinates_Reader: read: size(coordinates) is not"&
                //" correct.")
        end if

        do i_box = 1, size(this%components_coordinates, 2)
            call check_file_exists(coordinates(i_box)%string)
            open(newunit=coordinates_unit, recl=max_line_length, file=coordinates(i_box)%string, &
                status="old", action="read")
            read(coordinates_unit, *) comment_character, field, box_size
            call this%periodic_boxes(i_box)%set(box_size)
            call this%boxes_size_checker(i_box)%check()
            read(coordinates_unit, *) comment_character, field, nums_particles
            read(coordinates_unit, *) comment_character !components coordinates legend
            do i_component = 1, size(this%components_coordinates, 1)
                call this%components_coordinates(i_component, i_box)%reader%read(coordinates_unit, &
                     nums_particles(i_component))
            end do
            close(coordinates_unit)
        end do
    end subroutine Abstract_read

end module classes_complete_coordinates_reader
