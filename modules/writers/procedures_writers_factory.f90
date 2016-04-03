module procedures_writers_factory

use json_module, only: json_file
use data_wrappers_prefix, only: environment_prefix
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_number_to_string, only: Concrete_Number_to_String
use class_walls_potential, only: Abstract_Walls_Potential
use class_component_coordinates, only: Abstract_Component_Coordinates
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: mixture_set
use types_short_interactions_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Wrapper
use types_dipolar_interactions_wrapper, only: DES_Real_Pairs_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use class_line_writer, only: Abstract_Line_Writer, Concrete_Line_Writer, Null_Line_Writer
use class_triangle_writer, only: Concrete_Line_Selector, &
    Abstract_Triangle_Writer, Concrete_Triangle_Writer, Null_Triangle_Writer
use class_energy_writer, only: Abstract_Energy_Writer, Concrete_Energy_Writer, Null_Energy_Writer
use class_changes_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use class_component_coordinates_writer, only: Concrete_Coordinates_Writer_Selector, &
    Abstract_Coordinates_Writer, Concrete_Coordinates_Writer, Null_Coordinates_Writer
use types_writers_wrapper, only: Component_Writers_Wrapper, Writers_Wrapper
use procedures_property_inquirers, only: use_walls, component_exists, component_has_positions, &
    component_has_orientations, component_can_move, component_can_rotate, component_can_exchange, &
    component_is_dipolar, components_interact, component_interacts_with_wall

implicit none

private
public :: writers_create, writers_destroy

interface writers_create
    module procedure :: create_all
    module procedure :: create_components
    module procedure :: create_components_coordinates
    module procedure :: create_coordinates
    module procedure :: create_components_changes
    module procedure :: create_changes
    module procedure :: create_field
    module procedure :: create_walls
    module procedure :: create_switches
    module procedure :: create_short_energies
    module procedure :: create_dipolar_energies
    module procedure :: create_dipolar_mixture_energy
end interface writers_create

interface writers_destroy
    module procedure :: destroy_dipolar_mixture_energy
    module procedure :: destroy_triangle
    module procedure :: destroy_line
    module procedure :: destroy_changes
    module procedure :: destroy_coordinates
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface writers_destroy

contains

    subroutine create_all(writers, components, wall_pairs, short_pairs, dipolar_pairs, changes, &
        input_data, prefix)
        type(Writers_Wrapper), intent(out) :: writers
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        type(Pair_Potentials_Wrapper), intent(in) :: short_pairs(:)
        type(DES_Real_Pairs_Wrapper), intent(in) :: dipolar_pairs(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(components))

        call writers_create(writers%components, components, changes, input_data, prefix)
        call writers_create(writers%field, components, "field_energies.out")
        call writers_create(writers%walls, wall_pairs, "walls_energies.out")
        call writers_create(writers%switches, components, "switches.out")
        call writers_create(writers%short_energies, short_pairs, "short_energies.out")
        call writers_create(writers%dipolar_energies, dipolar_pairs, "dipolar_energies.out")
        call mixture_set(are_dipolar, components)
        call writers_create(writers%dipolar_mixture_energy, any(are_dipolar), &
            "dipolar_mixture_energy.out")
    end subroutine create_all

    subroutine destroy_all(writers)
        type(Writers_Wrapper), intent(inout) :: writers

        call writers_destroy(writers%dipolar_mixture_energy)
        call writers_destroy(writers%dipolar_energies)
        call writers_destroy(writers%short_energies)
        call writers_destroy(writers%switches)
        call writers_destroy(writers%walls)
        call writers_destroy(writers%field)
        call writers_destroy(writers%components)
    end subroutine destroy_all

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

    subroutine create_field(field, components, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: field
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        logical :: selector(size(components))
        integer :: i_component

        do i_component = 1, size(selector)
            selector(i_component) = component_is_dipolar(components(i_component)%dipolar_moments)
        end do

        if (any(selector)) then
            allocate(Concrete_Line_Writer :: field)
        else
            allocate(Null_Line_Writer :: field)
        end if
        call field%construct(filename, selector)
    end subroutine create_field

    subroutine create_walls(walls, wall_pairs, filename)
        class(Abstract_Line_Writer), allocatable, intent(out) :: walls
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        character(len=*), intent(in) :: filename

        logical :: selector(size(wall_pairs))
        integer :: i_component

        do i_component = 1, size(wall_pairs)
            selector(i_component) = component_interacts_with_wall(wall_pairs(i_component)%potential)
        end do

        if (any(selector)) then
            allocate(Concrete_Line_Writer :: walls)
        else
            allocate(Null_Line_Writer :: walls)
        end if
        call walls%construct(filename, selector)
    end subroutine create_walls

    subroutine destroy_line(line)
        class(Abstract_Line_Writer), allocatable, intent(inout) :: line

        if (allocated(line)) then
            call line%destroy()
            deallocate(line)
        end if
    end subroutine destroy_line

    subroutine create_switches(switches, components, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: switches
        type(Component_Wrapper), intent(in) :: components(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Line_Selector) :: selectors(size(components))
        logical :: some_components_exist, exist_ij
        integer :: i_component, j_component

        some_components_exist = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                exist_ij = component_exists(components(i_component)%number) .and. &
                    component_exists(components(j_component)%number)
                some_components_exist = some_components_exist .or. exist_ij
                selectors(j_component)%line(i_component) = exist_ij .and. i_component /= j_component
            end do
        end do
        if (some_components_exist) then
            allocate(Concrete_Triangle_Writer :: switches)
        else
            allocate(Null_Triangle_Writer :: switches)
        end if
        call switches%construct(filename, selectors)
        call deallocate_selectors(selectors)
    end subroutine create_switches

    subroutine create_short_energies(energies, pairs, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        type(Pair_Potentials_Wrapper), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Line_Selector) :: selectors(size(pairs))
        logical :: some_components_interact, interact_ij
        integer :: i_component, j_component

        some_components_interact = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                interact_ij = components_interact(pairs(j_component)%line(i_component)%potential)
                some_components_interact = some_components_interact .or. interact_ij
                selectors(j_component)%line(i_component) = interact_ij
            end do
        end do

        if (some_components_interact) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(filename, selectors)
        call deallocate_selectors(selectors)
    end subroutine create_short_energies

    subroutine create_dipolar_energies(energies, pairs, filename)
        class(Abstract_Triangle_Writer), allocatable, intent(out) :: energies
        type(DES_Real_Pairs_Wrapper), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Line_Selector) :: selectors(size(pairs))
        logical :: some_components_interact, interact_ij
        integer :: i_component, j_component

        some_components_interact = .false.
        do j_component = 1, size(selectors)
            allocate(selectors(j_component)%line(j_component))
            do i_component = 1, size(selectors(j_component)%line)
                interact_ij = components_interact(pairs(j_component)%line(i_component)%potential)
                some_components_interact = some_components_interact .or. interact_ij
                selectors(j_component)%line(i_component) = interact_ij
            end do
        end do

        if (some_components_interact) then
            allocate(Concrete_Triangle_Writer :: energies)
        else
            allocate(Null_Triangle_Writer :: energies)
        end if
        call energies%construct(filename, selectors)
        call deallocate_selectors(selectors)
    end subroutine create_dipolar_energies

    subroutine deallocate_selectors(selectors)
        type(Concrete_Line_Selector), intent(inout) :: selectors(:)

        integer :: i_component
        do i_component = size(selectors), 1, -1
            if (allocated(selectors(i_component)%line)) then
                deallocate(selectors(i_component)%line)
            end if
        end do
    end subroutine deallocate_selectors

    subroutine destroy_triangle(triangle)
        class(Abstract_Triangle_Writer), allocatable, intent(inout) :: triangle

        if (allocated(triangle)) then
            call triangle%destroy()
            deallocate(triangle)
        end if
    end subroutine destroy_triangle

    subroutine create_dipolar_mixture_energy(dipolar_mixture_energy, dipoles_exist, filename)
        class(Abstract_Energy_Writer), allocatable, intent(out) :: dipolar_mixture_energy
        logical, intent(in) :: dipoles_exist
        character(len=*), intent(in) :: filename

        if (dipoles_exist) then
            allocate(Concrete_Energy_Writer :: dipolar_mixture_energy)
        else
            allocate(Null_Energy_Writer :: dipolar_mixture_energy)
        end if
        call dipolar_mixture_energy%construct(filename)
    end subroutine create_dipolar_mixture_energy

    subroutine destroy_dipolar_mixture_energy(dipolar_mixture_energy)
        class(Abstract_Energy_Writer), allocatable, intent(inout) :: dipolar_mixture_energy

        if (allocated(dipolar_mixture_energy)) then
            call dipolar_mixture_energy%destroy()
            deallocate(dipolar_mixture_energy)
        end if
    end subroutine destroy_dipolar_mixture_energy

end module procedures_writers_factory
