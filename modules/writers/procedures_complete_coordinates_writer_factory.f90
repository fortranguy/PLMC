module procedures_complete_coordinates_writer_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use classes_component_coordinates_writer, only: Component_Coordinates_Writer_Selector
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper
use procedures_component_coordinates_writer_factory, only: component_coordinates_writer_create => &
    create, component_coordinates_writer_destroy => destroy
use classes_complete_coordinates_writer, only: Abstract_Complete_Coordinates_Writer, &
    Concrete_Complete_Coordinates_Writer, Null_Complete_Coordinates_Writer
use procedures_property_inquirers, only:  component_has_positions, component_has_orientations

implicit none

private
public :: create, destroy

contains

    subroutine create(coordinates, periodic_box, components, basename, generating_data, prefix)
        class(Abstract_Complete_Coordinates_Writer), allocatable, intent(out) :: coordinates
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: basename
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: period
        logical :: write_coordinates
        character(len=:), allocatable :: data_field
        logical :: data_found

        type(Component_Coordinates_Writer_Wrapper) :: components_coordinates(size(components))
        type(Component_Coordinates_Writer_Selector) :: components_selectors(size(components))

        data_field = prefix//"Coordinates.write"
        call generating_data%get(data_field, write_coordinates, data_found)
        call check_data_found(data_field, data_found)
        if (write_coordinates) then
            allocate(Concrete_Complete_Coordinates_Writer :: coordinates)
        else
            allocate(Null_Complete_Coordinates_Writer :: coordinates)
        end if
        if (write_coordinates) then
            data_field = prefix//"Coordinates.period"
            call generating_data%get(data_field, period, data_found)
            call check_data_found(data_field, data_found)
        end if

        call create_components(components_coordinates, components_selectors, components, &
            write_coordinates)
        call coordinates%construct(periodic_box, components_coordinates, components_selectors, &
            basename, period)
        call destroy_components(components_coordinates)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Complete_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

    subroutine create_components(components_coordinates, components_selectors, components, &
        write_coordinates)
        type(Component_Coordinates_Writer_Wrapper), intent(inout) :: components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), intent(inout) :: components_selectors(:)
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: write_coordinates

        integer :: i_component

        do i_component = 1, size(components_coordinates)
            associate(positions_i => components(i_component)%positions, &
                orientations_i => components(i_component)%orientations)
                components_selectors(i_component)%write_positions = &
                    component_has_positions(positions_i)
                components_selectors(i_component)%write_orientations = &
                    component_has_orientations(orientations_i)
                call component_coordinates_writer_create(components_coordinates(i_component)%&
                    writer, i_component, positions_i, orientations_i, &
                    components_selectors(i_component), write_coordinates)
            end associate
        end do
    end subroutine create_components

    subroutine destroy_components(components_coordinates)
        type(Component_Coordinates_Writer_Wrapper), intent(inout) :: components_coordinates(:)

        integer :: i_component

        do i_component = size(components_coordinates), 1, -1
            call component_coordinates_writer_destroy(components_coordinates(i_component)%writer)
        end do
    end subroutine destroy_components

end module procedures_complete_coordinates_writer_factory
