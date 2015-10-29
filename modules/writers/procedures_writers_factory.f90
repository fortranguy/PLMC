module procedures_writers_factory

use json_module, only: json_file
use data_wrappers_prefix, only: environment_prefix
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_number_to_string, only: Concrete_Number_to_String
use class_component_coordinates, only: Abstract_Component_Coordinates
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Pair_Potentials_Wrapper
use types_long_interactions_wrapper, only: Ewald_Real_Pairs_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use class_components_energes_writer, only: Concrete_Components_Energies_Selector, &
    Abstract_Components_Energies_Writer, Concrete_Components_Energies_Writer, &
    Null_Components_Energies_Writer
use class_changes_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use class_component_coordinates_writer, only: Concrete_Coordinates_Writer_Selector, &
    Abstract_Coordinates_Writer, Concrete_Coordinates_Writer, Null_Coordinates_Writer
use types_writers_wrapper, only: Component_Writers_Wrapper, Writers_Wrapper
use procedures_property_inquirers, only: use_walls, component_has_positions, &
    component_has_orientations, component_can_move, component_can_rotate, component_can_exchange, &
    component_is_dipolar, components_interact

implicit none

private
public :: writers_create, writers_destroy

interface writers_create
    module procedure :: create_all
    module procedure :: create_components_energies
    module procedure :: create_components
    module procedure :: create_components_changes
    module procedure :: create_changes
    module procedure :: create_components_coordinates
    module procedure :: create_coordinates
end interface writers_create

interface writers_destroy
    module procedure :: destroy_coordinates
    module procedure :: destroy_changes
    module procedure :: destroy_components
    module procedure :: destroy_components_energies
    module procedure :: destroy_all
end interface writers_destroy

contains

    subroutine create_all(writers, components, short_pairs, long_pairs, changes, input_data, &
        prefix)
        type(Writers_Wrapper), intent(out) :: writers
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Wrapper), intent(in) :: short_pairs(:)
        type(Ewald_Real_Pairs_Wrapper), intent(in) :: long_pairs(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        !todo: walls and field to add
        call writers_create(writers%short_energies, short_pairs, "short_energies.out")
        call writers_create(writers%long_energies, long_pairs, "long_energies.out")
        call writers_create(writers%components, components, changes, input_data, prefix)
    end subroutine create_all

    subroutine destroy_all(writers)
        type(Writers_Wrapper), intent(inout) :: writers

        call writers_destroy(writers%components)
        call writers_destroy(writers%long_energies)
        call writers_destroy(writers%short_energies)
    end subroutine destroy_all

    subroutine create_components_energies(energies, pairs, filename)
        class(Abstract_Components_Energies_Writer), allocatable, intent(out) :: energies
        class(*), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Components_Energies_Selector) :: selectors(size(pairs))
        logical :: interact, interact_ij
        integer :: j_component, i_component

        interact = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%with_components(j_component))
            do i_component = 1, size(selectors(j_component)%with_components)
                select type (pairs) ! several checks: not optimal
                    type is (Pair_Potentials_Wrapper)
                        interact_ij = components_interact(pairs(j_component)%&
                            with_components(i_component)%pair_potential)
                    type is (Ewald_Real_Pairs_Wrapper)
                        interact_ij = components_interact(pairs(j_component)%&
                            with_components(i_component)%real_pair)
                    class default
                        call error_exit("writers_create: create_components_energies: "//&
                            "pairs type unknown.")
                end select
                interact = interact .or. interact_ij
                selectors(j_component)%with_components(i_component) = interact_ij
            end do
        end do
        if (interact) then
            allocate(Concrete_Components_Energies_Writer :: energies)
        else
            allocate(Null_Components_Energies_Writer :: energies)
        end if
        call energies%construct(filename, selectors)
        call deallocate_selectors(selectors)
    end subroutine create_components_energies

    subroutine deallocate_selectors(selectors)
        type(Concrete_Components_Energies_Selector), intent(inout) :: selectors(:)

        integer :: i_component
        do i_component = size(selectors), 1, -1
            deallocate(selectors(i_component)%with_components)
        end do
    end subroutine deallocate_selectors

    subroutine destroy_components_energies(energies)
        class(Abstract_Components_Energies_Writer), allocatable, intent(inout) :: energies

        if (allocated(energies)) then
            call energies%destroy()
            deallocate(energies)
        end if
    end subroutine destroy_components_energies

    subroutine create_components(components, mixture_components, changes, input_data, prefix)
        type(Component_Writers_Wrapper), allocatable, intent(out) :: components(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        allocate(components(size(mixture_components)))
        call writers_create(components, mixture_components, input_data, prefix)
        call writers_create(components, changes)
    end subroutine create_components

    subroutine destroy_components(components)
        type(Component_Writers_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call writers_destroy(components(i_component)%changes)
                call writers_destroy(components(i_component)%coordinates)
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
            deallocate(data_field)
        end if
        do i_component = 1, size(components)
            associate(positions_i => mixture_components(i_component)%positions, &
                orientations_i => mixture_components(i_component)%orientations)
                selector_i%write_positions = component_has_positions(positions_i)
                selector_i%write_orientations = component_has_orientations(orientations_i)
                call writers_create(components(i_component)%coordinates, "coordinates_"//&
                    string%get(i_component), positions_i, orientations_i, selector_i, &
                    write_coordinates)
            end associate
        end do
    end subroutine create_components_coordinates

    subroutine create_coordinates(coordinates, basename, positions, orientations, &
        selector, write_coordinates)
        class(Abstract_Coordinates_Writer), allocatable, intent(out) :: coordinates
        character(len=*), intent(in) :: basename
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: selector
        logical, intent(in) :: write_coordinates

        if (write_coordinates) then
            allocate(Concrete_Coordinates_Writer :: coordinates)
        else
            allocate(Null_Coordinates_Writer :: coordinates)
        end if
        call coordinates%construct(basename, positions, orientations, selector)
    end subroutine create_coordinates

    subroutine destroy_coordinates(coordinates)
        class(Abstract_Coordinates_Writer), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_coordinates

    subroutine create_components_changes(components, changes)
        type(Component_Writers_Wrapper), intent(inout) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        type(Concrete_Changes_Selector) :: selector_i
        type(Concrete_Number_to_String) :: string

        integer :: i_component
        do i_component = 1, size(components)
            selector_i%write_positions = component_can_move(changes(i_component)%moved_positions)
            selector_i%write_rotations = component_can_rotate(changes(i_component)%&
                rotated_orientations)
            selector_i%write_exchanges = component_can_exchange(changes(i_component)%exchange)
            call writers_create(components(i_component)%changes, selector_i, "changes_"//&
                string%get(i_component)//".out")
        end do
    end subroutine create_components_changes

    subroutine create_changes(changes, selector, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes
        type(Concrete_Changes_Selector), intent(in) :: selector
        character(len=*), intent(in) :: filename

        if (selector%write_positions) then
            allocate(Concrete_Changes_Success_Writer :: changes)
        else
            allocate(Null_Changes_Success_Writer :: changes)
        end if
        call changes%construct(filename, selector)
    end subroutine create_changes

    subroutine destroy_changes(changes)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes

        if (allocated(changes)) then
            call changes%destroy()
            deallocate(changes)
        end if
    end subroutine destroy_changes

end module procedures_writers_factory

