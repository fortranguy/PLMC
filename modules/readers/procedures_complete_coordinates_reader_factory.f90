module procedures_complete_coordinates_reader_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use types_component_wrapper, only: Component_Wrapper
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper
use procedures_component_coordinates_reader_factory, only: &
    component_coordinates_reader_create => create, component_coordinates_reader_destroy => destroy
use classes_complete_coordinates_reader, only: Abstract_Complete_Coordinates_Reader, &
    Concrete_Complete_Coordinates_Reader

implicit none

private
public :: create, destroy

contains

    subroutine create(coordinates, periodic_box, box_size_checker, components)
        class(Abstract_Complete_Coordinates_Reader), allocatable, intent(out) :: coordinates
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Box_Size_Checker), intent(in) :: box_size_checker
        type(Component_Wrapper), intent(in) :: components(:)

        type(Component_Coordinates_Reader_wrapper), allocatable :: components_coordinates(:)

        allocate(Concrete_Complete_Coordinates_Reader :: coordinates)

        call component_coordinates_reader_create(components_coordinates, components)
        call coordinates%construct(periodic_box, box_size_checker, components_coordinates)
        call component_coordinates_reader_destroy(components_coordinates)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Complete_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_complete_coordinates_reader_factory
