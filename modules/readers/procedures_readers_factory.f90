module procedures_readers_factory

use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use class_component_number, only: Abstract_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates
use types_component_wrapper, only: Component_Wrapper
use class_box_size_reader, only: Abstract_Box_Size_Reader, Concrete_Box_Size_Reader
use class_component_coordinates_reader, only: Concrete_Coordinates_Reader_Selector, &
    Abstract_Coordinates_Reader, Concrete_Coordinates_Reader, Null_Coordinates_Reader
use types_readers_wrapper, only: Component_Readers_wrapper, Readers_Wrapper
use procedures_property_inquirers, only: component_has_positions, component_has_orientations

implicit none

private
public :: readers_create

interface readers_create
    module procedure :: create_all
    module procedure :: create_box_size
    module procedure :: create_components_coordinates
    module procedure :: create_coordinates
end interface readers_create

interface readers_destroy
    module procedure :: destroy_coordinates
    module procedure :: destroy_components_coordinates
    module procedure :: destroy_box_size
    module procedure :: destroy_all
end interface readers_destroy

contains

    subroutine create_all(readers, environment, components)
        type(Readers_Wrapper), intent(out) :: readers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)

        call readers_create(readers%box_size, environment%periodic_box)
        call readers_create(readers%components, components)
    end subroutine create_all

    subroutine destroy_all(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call readers_destroy(readers%components)
        call readers_destroy(readers%box_size)
    end subroutine destroy_all

    subroutine create_box_size(box_size, periodic_box)
        class(Abstract_Box_Size_Reader), allocatable, intent(out) :: box_size
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        allocate(Concrete_Box_Size_Reader :: box_size)
        call box_size%construct(periodic_box)
    end subroutine create_box_size

    subroutine destroy_box_size(box_size)
        class(Abstract_Box_Size_Reader), allocatable, intent(inout) :: box_size

        if (allocated(box_size)) then
            call box_size%destroy()
            deallocate(box_size)
        end if
    end subroutine destroy_box_size

    subroutine create_components_coordinates(components, mixture_components)
        type(Component_Readers_wrapper), intent(out) :: components(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:)

        type(Concrete_Coordinates_Reader_Selector) :: selector_i
        integer :: i_component

        do i_component = 1, size(components)
            associate(number_i => mixture_components(i_component)%number, &
                positions_i => mixture_components(i_component)%positions, &
                orientations_i => mixture_components(i_component)%orientations)
                selector_i%read_positions = component_has_positions(positions_i)
                selector_i%read_orientations = component_has_orientations(orientations_i)
                call readers_create(components(i_component)%coordinates, number_i, positions_i, &
                    orientations_i, selector_i)
            end associate
        end do
    end subroutine create_components_coordinates

    subroutine destroy_components_coordinates(components)
        type(Component_Readers_wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call readers_destroy(components(i_component)%coordinates)
            end do
            deallocate(components)
        end if
    end subroutine destroy_components_coordinates

    subroutine create_coordinates(coordinates, number, positions, orientations, selector)
        class(Abstract_Coordinates_Reader), allocatable, intent(out) :: coordinates
        class(Abstract_Component_Number), intent(in) :: number
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Concrete_Coordinates_Reader_Selector), intent(in) :: selector

        if (selector%read_positions) then
            allocate(Concrete_Coordinates_Reader :: coordinates)
        else
            allocate(Null_Coordinates_Reader :: coordinates)
        end if
        call coordinates%construct(number, positions, orientations, selector%read_orientations)
    end subroutine create_coordinates

    subroutine destroy_coordinates(coordinates)
        class(Abstract_Coordinates_Reader), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_coordinates

end module procedures_readers_factory
