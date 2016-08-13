module procedures_generating_writers_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use procedures_complete_coordinates_writer_factory, only: complete_coordinates_writer_create => &
    create, complete_coordinates_writer_destroy => destroy
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_factory, only: set_can_exchange
use procedures_real_writer_factory, only: real_writer_create => create, &
    real_writer_destroy => destroy
use procedures_line_writer_factory, only: line_writer_create => create, &
    line_writer_destroy => destroy
use procedures_triangle_writer_factory, only: triangle_writer_create => create, &
    triangle_writer_destroy => destroy
use procedures_square_writer_factory, only: square_writer_create => create_transmutations, &
    square_writer_destroy => destroy
use classes_changes_success_writer, only: Concrete_Changes_Selector
use procedures_changes_success_writer_factory, only: changes_success_writer_create => create, &
    changes_success_writer_destroy => destroy
use types_changes_success_writer_wrapper, only: Changes_Success_Writer_Wrapper
use types_generating_writers_wrapper, only: Generating_Writers_Wrapper
use procedures_property_inquirers, only:  component_can_translate, component_can_rotate, &
    component_can_exchange

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_all
    module procedure :: create_change_components
end interface create

interface destroy
    module procedure :: destroy_change_components
    module procedure :: destroy_all
end interface destroy

contains

    subroutine create_all(writers, environment, wall_pairs, components, short_pairs, &
        changes_components, generating_data, prefix)
        type(Generating_Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Line), intent(in) :: short_pairs(:)
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        logical, dimension(size(components)) :: can_exchange, are_dipolar

        call line_writer_create(writers%field, environment%external_field, components, &
            "field_energies.out")
        call line_writer_create(writers%walls, wall_pairs, "walls_energies.out")
        call set_can_exchange(can_exchange, components)
        call line_writer_create(writers%num_particles, can_exchange, "num_particles.out")
        call complete_coordinates_writer_create(writers%complete_coordinates, environment%&
            periodic_box, components, "coordinates", generating_data, prefix)
        call create(writers%components_changes, changes_components, components)
        call triangle_writer_create(writers%short_energies, short_pairs, "short_energies.out")
        call set_are_dipolar(are_dipolar, components)
        call triangle_writer_create(writers%dipolar_energies, are_dipolar, "dipolar_energies.out")
        call real_writer_create(writers%dipolar_mixture_energy, any(are_dipolar), &
            "dipolar_mixture_energy.out")
        call triangle_writer_create(writers%switches, components, "switches.out")
        call square_writer_create(writers%transmutations, components, "transmutations.out")
    end subroutine create_all

    subroutine destroy_all(writers)
        type(Generating_Writers_Wrapper), intent(inout) :: writers

        call square_writer_destroy(writers%transmutations)
        call triangle_writer_destroy(writers%switches)
        call real_writer_destroy(writers%dipolar_mixture_energy)
        call triangle_writer_destroy(writers%dipolar_energies)
        call triangle_writer_destroy(writers%short_energies)
        call destroy(writers%components_changes)
        call complete_coordinates_writer_destroy(writers%complete_coordinates)
        call line_writer_destroy(writers%num_particles)
        call line_writer_destroy(writers%walls)
        call line_writer_destroy(writers%field)
    end subroutine destroy_all

    subroutine create_change_components(components_changes, changes_components, components)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(out) :: components_changes(:)
        type(Changes_Component_Wrapper), intent(in) :: changes_components(:)
        type(Component_Wrapper), intent(in) :: components(:)

        type(Concrete_Changes_Selector) :: selector_i
        type(Concrete_Number_to_String) :: string
        integer :: i_component

        allocate(components_changes(size(changes_components)))
        do i_component = 1, size(components_changes)
            selector_i%write_translations = &
                component_can_translate(changes_components(i_component)%translated_positions)
            selector_i%write_rotations = component_can_rotate(changes_components(i_component)%&
                rotated_orientations)
            selector_i%write_exchanges = component_can_exchange(components(i_component)%&
                chemical_potential)
            call changes_success_writer_create(components_changes(i_component)%writer, selector_i, &
                "changes_"//string%get(i_component)//".out")
        end do
    end subroutine create_change_components

     subroutine destroy_change_components(components_changes)
        type(Changes_Success_Writer_Wrapper), allocatable, intent(inout) :: components_changes(:)

        integer :: i_component

        if (allocated(components_changes)) then
            do i_component = size(components_changes), 1, -1
                call changes_success_writer_destroy(components_changes(i_component)%writer)
            end do
            deallocate(components_changes)
        end if
    end subroutine destroy_change_components

end module procedures_generating_writers_factory
