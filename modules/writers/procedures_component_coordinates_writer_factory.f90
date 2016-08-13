module procedures_component_coordinates_writer_factory

use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_coordinates_writer, only: Component_Coordinates_Writer_Selector, &
    Abstract_Component_Coordinates_Writer, Concrete_Component_Coordinates_Writer, &
    Null_Component_Coordinates_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(writer, positions, orientations, selector, write_coordinates, basename)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(out) :: writer
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Component_Coordinates_Writer_Selector), intent(in) :: selector
        logical, intent(in) :: write_coordinates
        character(len=*), intent(in) :: basename

        if (write_coordinates .and. (selector%write_positions.or.selector%write_orientations)) then
            allocate(Concrete_Component_Coordinates_Writer :: writer)
        else
            allocate(Null_Component_Coordinates_Writer :: writer)
        end if
        call writer%construct(positions, orientations, selector, basename)
    end subroutine create

    subroutine destroy(writer)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy

end module procedures_component_coordinates_writer_factory
