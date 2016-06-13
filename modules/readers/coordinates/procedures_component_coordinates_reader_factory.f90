module procedures_component_coordinates_reader_factory

use classes_component_number, only: Abstract_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_coordinates_reader, only: Concrete_Component_Coordinates_Reader_Selector, &
    Abstract_Component_Coordinates_Reader, Concrete_Component_Coordinates_Reader, &
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

        if (selector%read_positions .or. selector%read_orientations) then
            allocate(Concrete_Component_Coordinates_Reader :: coordinates)
        else
            allocate(Null_Component_Coordinates_Reader :: coordinates)
        end if
        call coordinates%construct(number, positions, orientations, selector)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_component_coordinates_reader_factory
