module procedures_component_coordinates_writer_factory

use types_component_wrapper, only: Component_Wrapper
use classes_component_coordinates, only: Abstract_Component_Coordinates
use procedures_mixture_inquirers, only:  component_has_positions, component_has_orientations
use types_component_coordinates_writer_selector, only: Component_Coordinates_Writer_Selector
use classes_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer, &
    Concrete_Component_Coordinates_Writer, Null_Component_Coordinates_Writer
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains


    subroutine create_line(components_coordinates, components_selectors, components, &
        write_coordinates)
        type(Component_Coordinates_Writer_Wrapper), allocatable, intent(out) :: &
            components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), allocatable, intent(out) :: &
            components_selectors(:)
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: write_coordinates

        integer :: i_component

        allocate(components_coordinates(size(components)))
        allocate(components_selectors(size(components_coordinates)))
        do i_component = 1, size(components_coordinates)
            associate(positions_i => components(i_component)%positions, &
                orientations_i => components(i_component)%orientations)
                components_selectors(i_component)%write_positions = &
                    component_has_positions(positions_i)
                components_selectors(i_component)%write_orientations = &
                    component_has_orientations(orientations_i)
                call create(components_coordinates(i_component)%writer, i_component, positions_i, &
                    orientations_i, components_selectors(i_component), write_coordinates)
            end associate
        end do
    end subroutine create_line

    subroutine destroy_line(components_coordinates)
        type(Component_Coordinates_Writer_Wrapper), allocatable, intent(inout) :: &
            components_coordinates(:)

        integer :: i_component

        if (allocated(components_coordinates)) then
            do i_component = size(components_coordinates), 1, -1
                call destroy(components_coordinates(i_component)%writer)
            end do
            deallocate(components_coordinates)
        end if
    end subroutine destroy_line

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
