module procedures_complete_coordinates_writer_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use types_component_coordinates_writer_selector, only: Component_Coordinates_Writer_Selector
use types_component_coordinates_writer_wrapper, only: Component_Coordinates_Writer_Wrapper
use procedures_component_coordinates_writer_factory, only: component_coordinates_writer_create => &
    create, component_coordinates_writer_destroy => destroy
use classes_complete_coordinates_writer, only: Abstract_Complete_Coordinates_Writer, &
    Concrete_Complete_Coordinates_Writer, Null_Complete_Coordinates_Writer
use procedures_writers_inquirers, only: property_write_coordinates => write_coordinates

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

        type(Component_Coordinates_Writer_Wrapper), allocatable :: components_coordinates(:)
        type(Component_Coordinates_Writer_Selector), allocatable :: components_selectors(:)
        type(Component_Coordinates_Writer_Selector) :: components_selector

        write_coordinates = property_write_coordinates(generating_data, prefix)
        if (write_coordinates) then
            allocate(Concrete_Complete_Coordinates_Writer :: coordinates)
        else
            allocate(Null_Complete_Coordinates_Writer :: coordinates)
        end if
        if (write_coordinates) then
            data_field = prefix//"Coordinates.period"
            call generating_data%get(data_field, period, data_found)
            call check_data_found(data_field, data_found)
        else
            period = 0
        end if

        call component_coordinates_writer_create(components_coordinates, components_selectors, &
            components, write_coordinates)
        components_selector%write_positions = any(components_selectors%write_positions)
        components_selector%write_orientations = any(components_selectors%write_orientations)
        call coordinates%construct(periodic_box, components_coordinates, components_selector, &
            basename, period)
        call component_coordinates_writer_destroy(components_coordinates)
    end subroutine create

    subroutine destroy(coordinates)
        class(Abstract_Complete_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_complete_coordinates_writer_factory
