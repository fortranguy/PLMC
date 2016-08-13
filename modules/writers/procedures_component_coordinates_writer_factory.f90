module procedures_component_coordinates_writer_factory

use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_coordinates_writer, only: Concrete_Component_Coordinates_Writer_Selector, &
    Abstract_Component_Coordinates_Writer, Concrete_Component_Coordinates_Writer, Null_Component_Coordinates_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(coordinates, positions, orientations, selector, write_coordinates, basename)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(out) :: coordinates
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Concrete_Component_Coordinates_Writer_Selector), intent(in) :: selector
        logical, intent(in) :: write_coordinates
        character(len=*), intent(in) :: basename

        if (write_coordinates .and. (selector%write_positions.or.selector%write_orientations)) then
            allocate(Concrete_Component_Coordinates_Writer :: coordinates)
        else
            allocate(Null_Component_Coordinates_Writer :: coordinates)
        end if
        call coordinates%construct(positions, orientations, selector, basename)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_component_coordinates_writer_factory
