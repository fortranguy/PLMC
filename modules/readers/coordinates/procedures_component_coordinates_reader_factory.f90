module procedures_component_coordinates_reader_factory

use procedures_errors, only: error_exit
use classes_component_number, only: Abstract_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_coordinates_reader, only: Concrete_Component_Coordinates_Reader_Selector, &
    Abstract_Component_Coordinates_Reader, Concrete_Component_Coordinates_Reader, &
    Concrete_Component_Positions_Reader, Concrete_Component_Orientations_Reader, &
    Null_Component_Coordinates_Reader

implicit none

private
public :: create, destroy

contains

    subroutine create(coordinates, number, positions, orientations, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        class(Abstract_Component_Number), intent(in) :: number
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Concrete_Component_Coordinates_Reader_Selector), intent(in) :: selector

        call allocate(coordinates, selector)
        call construct(coordinates, number, positions, orientations)
    end subroutine create

    subroutine allocate(coordinates, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        type(Concrete_Component_Coordinates_Reader_Selector), intent(in) :: selector

        if (selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Coordinates_Reader :: coordinates)
        else if (selector%read_positions .and. .not.selector%read_orientations) then
            allocate(Concrete_Component_Positions_Reader :: coordinates)
        else if (.not.selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Orientations_Reader :: coordinates)
        else
            allocate(Null_Component_Coordinates_Reader :: coordinates)
        end if
    end subroutine allocate

    subroutine construct(coordinates, number, positions, orientations)
        class(Abstract_Component_Coordinates_Reader), intent(inout) :: coordinates
        class(Abstract_Component_Number), intent(in) :: number
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations

        select type (coordinates)
            type is (Concrete_Component_Coordinates_Reader)
                call coordinates%construct(number, positions, orientations)
            type is (Concrete_Component_Positions_Reader)
                call coordinates%construct(number, positions)
            type is (Concrete_Component_Orientations_Reader)
                call coordinates%construct(number, orientations)
            type is (Null_Component_Coordinates_Reader)
                call coordinates%destroy()
            class default
                call error_exit("procedures_component_coordinates_reader_factory: construct: "//&
                    "coordinates type unknown.")
        end select
    end subroutine construct

    subroutine destroy(coordinates)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_component_coordinates_reader_factory
