module classes_complete_coordinates_reader

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use data_strings, only: max_line_length, max_word_length
use procedures_checks, only: check_file_exists
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper
use procedures_component_coordinates_reader_factory, only: component_coordinates_reader_destroy => &
    destroy

implicit none

private

    type, abstract, public :: Abstract_Complete_Coordinates_Reader
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Box_Size_Checker), pointer :: box_size_checker => null()
        type(Component_Coordinates_Reader_wrapper), allocatable :: components_coordinates(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: read => Abstract_read
    end type Abstract_Complete_Coordinates_Reader

    type, extends(Abstract_Complete_Coordinates_Reader), public :: &
        Concrete_Complete_Coordinates_Reader
    end type Concrete_Complete_Coordinates_Reader

contains

    subroutine Abstract_construct(this, periodic_box, box_size_checker, components_coordinates)
        class(Abstract_Complete_Coordinates_Reader), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Box_Size_Checker), target, intent(in) :: box_size_checker
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_coordinates(:)

        this%periodic_box => periodic_box
        this%box_size_checker => box_size_checker
        allocate(this%components_coordinates, source=components_coordinates)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Complete_Coordinates_Reader), intent(out) :: this

        call component_coordinates_reader_destroy(this%components_coordinates)
        this%box_size_checker => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_read(this, coordinates_filename)
        class(Abstract_Complete_Coordinates_Reader), intent(in) :: this
        character(len=*), intent(in) :: coordinates_filename

        real(DP) :: box_size(num_dimensions)
        integer :: nums_particles(size(this%components_coordinates)), i_component
        character(len=1) :: comment_character
        character(len=max_word_length) :: field
        integer :: coordinates_unit

        call check_file_exists(coordinates_filename)
        open(newunit=coordinates_unit, recl=max_line_length, file=coordinates_filename, &
            status="old", action="read")
        read(coordinates_unit, *) comment_character, field, box_size
        call this%periodic_box%set(box_size)
        call this%box_size_checker%check()
        read(coordinates_unit, *) comment_character, field, nums_particles
        read(coordinates_unit, *) comment_character !components coordinates legend
        do i_component = 1, size(this%components_coordinates)
            call this%components_coordinates(i_component)%reader%read(coordinates_unit, &
                 nums_particles(i_component))
        end do
        close(coordinates_unit)
    end subroutine Abstract_read

end module classes_complete_coordinates_reader
