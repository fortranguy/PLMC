module procedures_complete_coordinates_writer_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use types_string_wrapper, only: String_Wrapper
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only:  component_has_positions, component_has_orientations
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

    subroutine create(coordinates, paths, basename, periodic_boxes, components, generating_data, &
        prefix)
        class(Abstract_Complete_Coordinates_Writer), allocatable, intent(out) :: coordinates
        type(String_Wrapper), intent(in) :: paths(:)
        character(len=*), intent(in) :: basename
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(Component_Wrapper), intent(in) :: components(:, :)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: period
        logical :: write_coordinates
        character(len=:), allocatable :: data_field
        logical :: data_found

        type(Component_Coordinates_Writer_Wrapper) :: components_coordinates(size(components, 1), &
            size(components, 2))
        type(Component_Coordinates_Writer_Selector) :: coordinates_selector

        write_coordinates = property_write_coordinates(generating_data, prefix)
        if (write_coordinates) then
            data_field = prefix//"Coordinates.period"
            call generating_data%get(data_field, period, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Complete_Coordinates_Writer :: coordinates)
        else
            period = 0
            allocate(Null_Complete_Coordinates_Writer :: coordinates)
        end if

        call set_selector(coordinates_selector, components)
        call component_coordinates_writer_create(components_coordinates, components, &
            write_coordinates)
        call coordinates%construct(paths, basename, periodic_boxes, components_coordinates, &
            coordinates_selector, period)
        call component_coordinates_writer_destroy(components_coordinates)
    end subroutine create

    subroutine set_selector(coordinates_selector, components)
        type(Component_Coordinates_Writer_Selector) :: coordinates_selector
        type(Component_Wrapper), intent(in) :: components(:, :)

        logical, dimension(size(components, 1), size(components, 2)) :: have_positions, &
            have_orientations
        integer :: i_box, i_component

        do i_box = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                have_positions(i_component, i_box) = &
                    component_has_positions(components(i_component, i_box)%positions)
                have_orientations(i_component, i_box) = &
                    component_has_orientations(components(i_component, i_box)%orientations)
            end do
        end do
        coordinates_selector%write_positions = any(have_positions)
        coordinates_selector%write_orientations = any(have_orientations)
    end subroutine set_selector

    subroutine destroy(coordinates)
        class(Abstract_Complete_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_complete_coordinates_writer_factory
