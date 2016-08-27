module procedures_component_coordinates_reader_factory

use procedures_errors, only: error_exit
use classes_component_number, only: Abstract_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_has_positions, component_has_orientations
use types_component_coordinates_reader_selector, only: Component_Coordinates_Reader_Selector
use classes_component_coordinates_reader, only: Abstract_Component_Coordinates_Reader, &
    Concrete_Component_Coordinates_Reader, Concrete_Component_Positions_Reader, &
    Concrete_Component_Orientations_Reader, Null_Component_Coordinates_Reader
use types_component_coordinates_reader_wrapper, only: Component_Coordinates_Reader_wrapper

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

    subroutine create_line(components_coordinates, components)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(out) :: &
            components_coordinates(:)
        type(Component_Wrapper), intent(in) :: components(:)

        type(Component_Coordinates_Reader_Selector) :: selector_i
        integer :: i_component

        allocate(components_coordinates(size(components)))
        do i_component = 1, size(components_coordinates)
            associate(number_i => components(i_component)%number, &
                positions_i => components(i_component)%positions, &
                orientations_i => components(i_component)%orientations)
                selector_i%read_positions = component_has_positions(positions_i)
                selector_i%read_orientations = component_has_orientations(orientations_i)
                call create(components_coordinates(i_component)%reader, number_i, positions_i, &
                    orientations_i, selector_i)
            end associate
        end do
    end subroutine create_line

    subroutine destroy_line(components_coordinates)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(inout) :: &
            components_coordinates(:)

        integer :: i_component

        if (allocated(components_coordinates)) then
            do i_component = size(components_coordinates), 1, -1
                call destroy(components_coordinates(i_component)%reader)
            end do
            deallocate(components_coordinates)
        end if
    end subroutine destroy_line

    subroutine create_element(coordinates, number, positions, orientations, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        class(Abstract_Component_Number), intent(in) :: number
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        call allocate_element(coordinates, selector)
        call construct_element(coordinates, number, positions, orientations)
    end subroutine create_element

    subroutine allocate_element(coordinates, selector)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(out) :: coordinates
        type(Component_Coordinates_Reader_Selector), intent(in) :: selector

        if (selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Coordinates_Reader :: coordinates)
        else if (selector%read_positions .and. .not.selector%read_orientations) then
            allocate(Concrete_Component_Positions_Reader :: coordinates)
        else if (.not.selector%read_positions .and. selector%read_orientations) then
            allocate(Concrete_Component_Orientations_Reader :: coordinates)
        else
            allocate(Null_Component_Coordinates_Reader :: coordinates)
        end if
    end subroutine allocate_element

    subroutine construct_element(coordinates, number, positions, orientations)
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
    end subroutine construct_element

    subroutine destroy_element(coordinates)
        class(Abstract_Component_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_element

end module procedures_component_coordinates_reader_factory
