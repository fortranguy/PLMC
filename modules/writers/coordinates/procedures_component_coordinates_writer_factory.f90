module procedures_component_coordinates_writer_factory

use types_component_wrapper, only: Component_Wrapper
use classes_component_coordinates, only: Abstract_Component_Coordinates
use procedures_mixture_inquirers, only:  component_has_positions, component_has_orientations
use types_component_coordinates_writer_selector, only: Component_Coordinates_Writer_Selector
use classes_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer, &
    Concrete_Component_Coordinates_Writer, Null_Component_Coordinates_Writer, &
    Component_Coordinates_Writer_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_rectangle
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_rectangle
end interface destroy

contains


    subroutine create_rectangle(components_coordinates, components, write_coordinates)
        type(Component_Coordinates_Writer_Wrapper), intent(inout) :: components_coordinates(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)
        logical, intent(in) :: write_coordinates

        type(Component_Coordinates_Writer_Selector) :: selector_i
        integer :: i_box, i_component

        do i_box = 1, size(components_coordinates, 2)
            do i_component = 1, size(components_coordinates, 1)
                associate(positions_i => components(i_component, i_box)%positions, &
                    orientations_i => components(i_component, i_box)%orientations)
                    selector_i%write_positions = component_has_positions(positions_i)
                    selector_i%write_orientations = component_has_orientations(orientations_i)
                    call create(components_coordinates(i_component, i_box)%writer, i_component, &
                        positions_i, orientations_i, selector_i, write_coordinates)
                end associate
            end do
        end do
    end subroutine create_rectangle

    subroutine destroy_rectangle(components_coordinates)
        type(Component_Coordinates_Writer_Wrapper), intent(inout) :: components_coordinates(:, :)

        integer :: i_box, i_component

        do i_box = size(components_coordinates, 2), 1, -1
            do i_component = size(components_coordinates, 1), 1, -1
                call destroy(components_coordinates(i_component, i_box)%writer)
            end do
        end do
    end subroutine destroy_rectangle

    subroutine create_element(coordinates, i_component, positions, orientations, selector, &
        write_coordinates)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(out) :: coordinates
        integer, intent(in) :: i_component
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Component_Coordinates_Writer_Selector), intent(in) :: selector
        logical, intent(in) :: write_coordinates

        if (write_coordinates .and. (selector%write_positions.or.selector%write_orientations)) then
            allocate(Concrete_Component_Coordinates_Writer :: coordinates)
        else
            allocate(Null_Component_Coordinates_Writer :: coordinates)
        end if
        call coordinates%construct(i_component, positions, orientations, selector)
    end subroutine create_element

    subroutine destroy_element(coordinates)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_element

end module procedures_component_coordinates_writer_factory
