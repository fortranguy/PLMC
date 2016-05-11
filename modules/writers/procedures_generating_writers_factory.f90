module procedures_generating_writers_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use classes_changes_success_writer, only: Concrete_Changes_Selector
use procedures_changes_success_writer_factory, only: changes_success_writer_create => create, &
    changes_success_writer_destroy => destroy
use classes_component_coordinates_writer, only: Concrete_Coordinates_Writer_Selector
use procedures_component_coordinates_writer_factory, only: &
    component_coordinates_writer_create => create, component_coordinates_writer_destroy => destroy
use types_component_writers_wrapper, only: Component_Writers_Wrapper
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use procedures_property_inquirers, only:  component_has_positions, component_has_orientations, &
    component_can_move, component_can_rotate, component_can_exchange

implicit none

private
public :: generating_writers_create, generating_writers_destroy

interface generating_writers_create
    module procedure :: create_all
    module procedure :: create_components
    module procedure :: create_components_coordinates
    module procedure :: create_change_components
end interface generating_writers_create

interface generating_writers_destroy
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface generating_writers_destroy

contains

    subroutine create_all(writers, environment, wall_pairs, components, change_components, &
        short_pairs, input_data, prefix)
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(components))

        call generating_writers_create(writers%components, components, change_components, &
            input_data, prefix)
        call line_writer_create(writers%field, environment%external_field, components, &
            "field_energies.out")
        call line_writer_create(writers%walls, wall_pairs, "walls_energies.out")
        call triangle_writer_create(writers%switches, components, "switches.out")
        call triangle_writer_create(writers%short_energies, short_pairs, "short_energies.out")
        call set_are_dipolar(are_dipolar, components)
        call triangle_writer_create(writers%dipolar_energies, are_dipolar, "dipolar_energies.out")
        call real_writer_create(writers%dipolar_mixture_energy, any(are_dipolar), &
            "dipolar_mixture_energy.out")
    end subroutine create_all

    subroutine destroy_all(writers)
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call real_writer_destroy(writers%dipolar_mixture_energy)
        call triangle_writer_destroy(writers%dipolar_energies)
        call triangle_writer_destroy(writers%short_energies)
        call triangle_writer_destroy(writers%switches)
        call line_writer_destroy(writers%walls)
        call line_writer_destroy(writers%field)
        call generating_writers_destroy(writers%components)
    end subroutine destroy_all

    subroutine create_components(components, mixture_components, change_components, input_data, &
        prefix)
        type(Component_Writers_Wrapper), allocatable, intent(out) :: components(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        allocate(components(size(mixture_components)))
        call generating_writers_create(components, mixture_components, input_data, prefix)
        call generating_writers_create(components, change_components)
    end subroutine create_components

    subroutine destroy_components(components)
        type(Component_Writers_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call changes_success_writer_destroy(components(i_component)%changes)
                call component_coordinates_writer_destroy(components(i_component)%coordinates)
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine create_components_coordinates(components, mixture_components, input_data, prefix)
        type(Component_Writers_Wrapper), intent(inout) :: components(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, write_coordinates
        type(Concrete_Coordinates_Writer_Selector) :: selector_i
        integer :: i_component
        type(Concrete_Number_to_String) :: string

        data_field = prefix//"Coordinates.write"
        call input_data%get(data_field, write_coordinates, data_found)
        call check_data_found(data_field, data_found)

        if (write_coordinates) then
            data_field = prefix//"Coordinates.period"
            call input_data%get(data_field, selector_i%period, data_found)
            call check_data_found(data_field, data_found)
        end if
        do i_component = 1, size(components)
            associate(positions_i => mixture_components(i_component)%positions, &
                orientations_i => mixture_components(i_component)%orientations)
                selector_i%write_positions = component_has_positions(positions_i)
                selector_i%write_orientations = component_has_orientations(orientations_i)
                call component_coordinates_writer_create(components(i_component)%coordinates, &
                    positions_i, orientations_i, selector_i, write_coordinates, &
                    "coordinates_"//string%get(i_component))
            end associate
        end do
    end subroutine create_components_coordinates

    subroutine create_change_components(components, change_components)
        type(Component_Writers_Wrapper), intent(inout) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        type(Concrete_Changes_Selector) :: selector_i
        type(Concrete_Number_to_String) :: string

        integer :: i_component
        do i_component = 1, size(components)
            selector_i%write_positions = component_can_move(change_components(i_component)%&
                moved_positions)
            selector_i%write_rotations = component_can_rotate(change_components(i_component)%&
                rotated_orientations)
            selector_i%write_exchanges = component_can_exchange(change_components(i_component)%&
                exchange)
            call changes_success_writer_create(components(i_component)%changes, selector_i, &
                "changes_"//string%get(i_component)//".out")
        end do
    end subroutine create_change_components

end module procedures_generating_writers_factory
