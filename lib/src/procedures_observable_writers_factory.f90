module procedures_observable_writers_factory

use json_module, only: json_file
use data_wrappers_prefix, only: environment_prefix
use procedures_checks, only: check_data_found
use class_component_number, only: Abstract_Component_Number
use class_component_diameter, only: Abstract_Component_Diameter
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_component_exchange, only: Abstract_Component_Exchange
use class_moved_positions, only: Abstract_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_component_energy_writer, only: Concrete_Energy_Writer_Selector, &
    Abstract_Component_Energy_Writer, Concrete_Component_Energy_Writer, Null_Component_Energy_Writer
use class_inter_energy_writer, only: Concrete_Inter_Energy_Writer_Selector, &
    Abstract_Inter_Energy_Writer, Concrete_Inter_Energy_Writer, Null_Inter_Energy_Writer
use class_changes_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use class_component_coordinates_writer, only: Concrete_Coordinates_Writer_Selector, &
    Abstract_Component_Coordinates_Writer, Concrete_Component_Coordinates_Writer, &
    Null_Component_Coordinates_Writer
use procedures_property_inquirers, only: use_walls, component_has_positions, &
    component_has_orientations, component_can_move, component_can_rotate, component_can_exchange, &
    component_is_dipolar

implicit none

private
public :: observable_writers_factory_create, observable_writers_factory_destroy

interface observable_writers_factory_create
    module procedure :: allocate_and_construct_energy
    module procedure :: allocate_and_construct_inter_energy
    module procedure :: allocate_and_construct_changes
    module procedure :: allocate_and_construct_coordinates
end interface observable_writers_factory_create

interface observable_writers_factory_destroy
    module procedure :: destroy_and_deallocate_coordinates
    module procedure :: destroy_and_deallocate_changes
    module procedure :: destroy_and_deallocate_inter_energy
    module procedure :: destroy_and_deallocate_component_energy
end interface observable_writers_factory_destroy

contains

    subroutine allocate_and_construct_energy(component_energy_writer, use_walls, component_exists, &
        component_is_dipolar, filename)
        class(Abstract_Component_Energy_Writer), allocatable, intent(out) :: component_energy_writer
        logical, intent(in) :: use_walls, component_exists, component_is_dipolar
        character(len=*), intent(in) :: filename

        type(Concrete_Energy_Writer_Selector) :: energy_selector

        if (component_exists) then
            allocate(Concrete_Component_Energy_Writer :: component_energy_writer)
        else
            allocate(Null_Component_Energy_Writer :: component_energy_writer)
        end if
        energy_selector%write_walls = use_walls .and. component_exists
        energy_selector%write_long = component_is_dipolar
        call component_energy_writer%construct(filename, energy_selector)
    end subroutine allocate_and_construct_energy

    subroutine allocate_and_construct_inter_energy(inter_energy_writer, mixture_exists, &
        mixture_is_dipolar, filename)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(out) :: inter_energy_writer
        logical, intent(in) :: mixture_exists, mixture_is_dipolar
        character(len=*), intent(in) :: filename

        type(Concrete_Inter_Energy_Writer_Selector) :: energy_selector

        if (mixture_exists) then
            allocate(Concrete_Inter_Energy_Writer :: inter_energy_writer)
        else
            allocate(Null_Inter_Energy_Writer :: inter_energy_writer)
        end if
        energy_selector%write_long = mixture_is_dipolar
        call inter_energy_writer%construct(filename, energy_selector)
    end subroutine allocate_and_construct_inter_energy

    subroutine allocate_and_construct_changes(changes_success_writer, moved_positions, &
        rotated_orientations, component_exchange, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes_success_writer
        class(Abstract_Moved_Positions), intent(in) :: moved_positions
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations
        class(Abstract_Component_Exchange), intent(in) :: component_exchange
        character(len=*), intent(in) :: filename

        type(Concrete_Changes_Selector) :: changes_selector

        if (component_can_move(moved_positions)) then
            allocate(Concrete_Changes_Success_Writer :: changes_success_writer)
        else
            allocate(Null_Changes_Success_Writer :: changes_success_writer)
        end if
        changes_selector%write_rotations = component_can_rotate(rotated_orientations)
        changes_selector%write_exchanges = component_can_exchange(component_exchange)
        call changes_success_writer%construct(filename, changes_selector)
    end subroutine allocate_and_construct_changes

    subroutine allocate_and_construct_coordinates(component_coordinates_writer, basename, &
        component_positions, component_orientations, input_data, prefix)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(out) :: &
            component_coordinates_writer
        character(len=*), intent(in) :: basename
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, write_coordinates
        type(Concrete_Coordinates_Writer_Selector) :: coordinates_selector

        data_field = prefix//"Coordinates.write"
        call input_data%get(data_field, write_coordinates, data_found)
        call check_data_found(data_field, data_found)

        if (write_coordinates .and. component_has_positions(component_positions)) then
            data_field = prefix//"Coordinates.period"
            call input_data%get(data_field, coordinates_selector%period, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Coordinates_Writer :: component_coordinates_writer)
        else
            allocate(Null_Component_Coordinates_Writer :: component_coordinates_writer)
        end if
        coordinates_selector%write_orientations = &
            component_has_orientations(component_orientations)
        call component_coordinates_writer%construct(basename, component_positions, &
            component_orientations, coordinates_selector)
        deallocate(data_field)
    end subroutine allocate_and_construct_coordinates

    subroutine destroy_and_deallocate_component_energy(component_energy_writer)
        class(Abstract_Component_Energy_Writer), allocatable, intent(inout) :: &
            component_energy_writer

        call component_energy_writer%destroy()
        if (allocated(component_energy_writer)) deallocate(component_energy_writer)
    end subroutine destroy_and_deallocate_component_energy

    subroutine destroy_and_deallocate_inter_energy(inter_energy_writer)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(inout) :: inter_energy_writer

        call inter_energy_writer%destroy()
        if (allocated(inter_energy_writer)) deallocate(inter_energy_writer)
    end subroutine destroy_and_deallocate_inter_energy

    subroutine destroy_and_deallocate_changes(changes_success_writer)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes_success_writer

        call changes_success_writer%destroy()
        if (allocated(changes_success_writer)) deallocate(changes_success_writer)
    end subroutine destroy_and_deallocate_changes

    subroutine destroy_and_deallocate_coordinates(component_coordinates_writer)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(inout) :: &
            component_coordinates_writer

        call component_coordinates_writer%destroy()
        if (allocated(component_coordinates_writer)) deallocate(component_coordinates_writer)
    end subroutine destroy_and_deallocate_coordinates

end module procedures_observable_writers_factory

