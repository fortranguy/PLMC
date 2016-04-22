module procedures_coordinates_reader_factory

use types_readers_wrapper, only: Component_Coordinates_Reader_wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_component_coordinates_reader, only: Concrete_Component_Coordinates_Reader_Selector
use procedures_component_coordinates_reader_factory, only: &
    component_coordinates_reader_create => create, component_coordinates_reader_destroy => destroy
use procedures_property_inquirers, only: component_has_positions, component_has_orientations

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_components
    module procedure :: component_coordinates_reader_create
end interface create

interface destroy
    module procedure :: component_coordinates_reader_destroy
    module procedure :: destroy_components
end interface destroy

contains

    subroutine create_components(components, mixture_components)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(out) :: components(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:)

        type(Concrete_Component_Coordinates_Reader_Selector) :: selector_i
        integer :: i_component

        allocate(components(size(mixture_components)))
        do i_component = 1, size(components)
            associate(number_i => mixture_components(i_component)%number, &
                positions_i => mixture_components(i_component)%positions, &
                orientations_i => mixture_components(i_component)%orientations)
                selector_i%read_positions = component_has_positions(positions_i)
                selector_i%read_orientations = component_has_orientations(orientations_i)
                call create(components(i_component)%coordinates, number_i, positions_i, &
                    orientations_i, selector_i)
            end associate
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(Component_Coordinates_Reader_wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call destroy(components(i_component)%coordinates)
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

end module procedures_coordinates_reader_factory
