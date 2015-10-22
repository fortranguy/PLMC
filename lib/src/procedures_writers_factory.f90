module procedures_writers_factory

use json_module, only: json_file
use data_wrappers_prefix, only: environment_prefix
use procedures_checks, only: check_data_found
use class_component_number, only: Abstract_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_short_potentials_wrapper, only: Pair_Potentials_Wrapper
use types_ewalds_wrapper, only: Ewald_Real_Pairs_Wrapper
use class_component_exchange, only: Abstract_Component_Exchange
use class_changed_coordinates, only: Abstract_Changed_Coordinates
use types_observables_wrapper, only: Observables_Wrapper
use class_inter_energes_writer, only: Concrete_Inter_Energy_Selector, Abstract_Inter_Energy_Writer
use class_changes_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer
use class_component_coordinates_writer, only: Concrete_Coordinates_Writer_Selector, &
    Abstract_Component_Coordinates_Writer, Concrete_Component_Coordinates_Writer, &
    Null_Component_Coordinates_Writer
use types_writers_wrapper, only: Writers_Wrapper
use procedures_property_inquirers, only: use_walls, component_has_positions, &
    component_has_orientations, component_can_move, component_can_rotate, component_can_exchange, &
    component_is_dipolar, components_interact

implicit none

private
public :: writers_create, writers_destroy

interface writers_create
    module procedure :: create_all
    module procedure :: create_short_inter
    module procedure :: create_long_inter
    module procedure :: create_changes
    module procedure :: create_component_coordinates
end interface writers_create

interface writers_destroy
    module procedure :: destroy_component_coordinates
    module procedure :: destroy_changes
    module procedure :: destroy_inter
    module procedure :: destroy_all
end interface writers_destroy

contains

    subroutine create_all(writers, components, short_pairs, ewald_pairs, observables, input_data, &
        prefix)
        type(Writers_Wrapper), intent(out) :: writers
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Wrapper), intent(in) :: short_pairs(:)
        type(Ewald_Real_Pairs_Wrapper), intent(in) :: ewald_pairs(:)
        type(Observables_Wrapper), intent(in) :: observables
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call writers_create(writers%short_inter, short_pairs, "short_potentials.out")
        call writers_create(writers%long_inter, ewald_pairs, "long_potentials.out")
    end subroutine create_all

    subroutine destroy_all(writers)
        type(Writers_Wrapper), intent(inout) :: writers

        call writers_destroy(writers%short_inter)
    end subroutine destroy_all

    subroutine create_short_inter(inter, pairs, filename)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(out) :: inter
        type(Pair_Potentials_Wrapper), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Inter_Energy_Selector) :: selector(size(pairs))
        integer :: j_component, i_component

        do j_component = 1, size(selector)
            allocate(selector(j_component)%with_components(j_component))
            do i_component = 1, size(selector(j_component)%with_components)
                associate(interact_ij => pairs(j_component)%with_components(i_component)%&
                    pair_potential)
                    selector(j_component)%with_components(i_component) = &
                        components_interact(interact_ij)
                end associate
            end do
        end do
        call inter%construct(filename, selector)
        call deallocate_selector(selector)
    end subroutine create_short_inter

    subroutine create_long_inter(inter, pairs, filename)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(out) :: inter
        type(Ewald_Real_Pairs_Wrapper), intent(in) :: pairs(:)
        character(len=*), intent(in) :: filename

        type(Concrete_Inter_Energy_Selector) :: selector(size(pairs))
        integer :: j_component, i_component

        do j_component = 1, size(selector)
            allocate(selector(j_component)%with_components(j_component))
            do i_component = 1, size(selector(j_component)%with_components)
                associate(interact_ij => pairs(j_component)%with_components(i_component)%&
                    real_pair)
                    selector(j_component)%with_components(i_component) = &
                        components_interact(interact_ij)
                end associate
            end do
        end do
        call inter%construct(filename, selector)
        call deallocate_selector(selector)
    end subroutine create_long_inter

    subroutine deallocate_selector(selector)
        type(Concrete_Inter_Energy_Selector), intent(inout) :: selector(:)

        integer :: i_component
        do i_component = size(selector), 1, -1
            deallocate(selector(i_component)%with_components)
        end do
    end subroutine deallocate_selector

    subroutine destroy_inter(inter)
        class(Abstract_Inter_Energy_Writer), allocatable, intent(inout) :: inter

        if (allocated(inter)) then
            call inter%destroy()
            deallocate(inter)
        end if
    end subroutine destroy_inter

    subroutine create_changes(changes_success_writer, changes_selector, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes_success_writer
        type(Concrete_Changes_Selector), intent(in) :: changes_selector
        character(len=*), intent(in) :: filename

        if (changes_selector%write_positions) then
            allocate(Concrete_Changes_Success_Writer :: changes_success_writer)
        else
            allocate(Null_Changes_Success_Writer :: changes_success_writer)
        end if
        call changes_success_writer%construct(filename, changes_selector)
    end subroutine create_changes

    subroutine destroy_changes(changes_success_writer)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes_success_writer

        if (allocated(changes_success_writer)) then
            call changes_success_writer%destroy()
            deallocate(changes_success_writer)
        end if
    end subroutine destroy_changes

    subroutine create_coordinates(input_data, prefix)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, write_coordinates
        type(Concrete_Coordinates_Writer_Selector) :: selector

        data_field = prefix//"Coordinates.write"
        call input_data%get(data_field, write_coordinates, data_found)
        call check_data_found(data_field, data_found)

        data_field = prefix//"Coordinates.period"
        call input_data%get(data_field, selector%period, data_found)
        call check_data_found(data_field, data_found)
        deallocate(data_field)
        !to do
    end subroutine create_coordinates

    subroutine create_component_coordinates(writer, basename, positions, orientations, selector)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(out) :: writer
        character(len=*), intent(in) :: basename
        class(Abstract_Component_Coordinates), intent(in) :: positions, orientations
        type(Concrete_Coordinates_Writer_Selector), intent(in) :: selector

        if (selector%write_positions) then
            allocate(Concrete_Component_Coordinates_Writer :: writer)
        else
            allocate(Null_Component_Coordinates_Writer :: writer)
        end if
        call writer%construct(basename, positions, orientations, selector)
    end subroutine create_component_coordinates

    subroutine destroy_component_coordinates(writer)
        class(Abstract_Component_Coordinates_Writer), allocatable, intent(inout) ::writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy_component_coordinates

end module procedures_writers_factory

