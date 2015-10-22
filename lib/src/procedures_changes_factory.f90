module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_plmc_iterations, only: num_tuning_steps
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_changed_coordinates, only: Abstract_Changed_Coordinates, Null_Changed_Coordinates
use class_moved_positions, only: Concrete_Moved_Positions
use class_rotated_orientations, only: Concrete_Rotated_Orientations
use class_change_tuner, only: Concrete_Change_Tuner_Parameters, Abstract_Change_Tuner, &
    Concrete_Change_Tuner, Null_Change_Tuner
use class_component_exchange, only: Abstract_Component_Exchange, Concrete_Component_Exchange, &
    Null_Component_Exchange
use module_change_tuning, only: Concrete_Tuning_Parameters
use types_component_wrapper, only: Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_property_inquirers, only: component_has_positions, component_has_orientations, &
    component_can_move, component_can_rotate, component_can_exchange

implicit none

private
public :: changes_create, changes_destroy, changes_create_move_tuner, changes_create_rotation_tuner

interface changes_create
    module procedure :: create_all
    module procedure :: create_moved_positions
    module procedure :: create_rotated_orientations
    module procedure :: create_component_exchange
end interface changes_create

interface changes_destroy
    module procedure :: destroy_component_exchange
    module procedure :: destroy_change_tuner
    module procedure :: destroy_changed_coordinates
    module procedure :: destroy_all
end interface changes_destroy

contains

    subroutine create_all(changes, periodic_box, component, input_data, prefix)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: component
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call changes_create(changes%moved_positions, component%positions, periodic_box, &
            input_data, prefix)
        call changes_create_move_tuner(changes%move_tuner, changes%moved_positions, &
            input_data, prefix)
        call changes_create(changes%rotated_orientations, component%orientations, &
            input_data, prefix)
        call changes_create_rotation_tuner(changes%rotation_tuner, &
            changes%rotated_orientations, input_data, prefix)
        call changes_create(changes%component_exchange, component)
    end subroutine create_all

    subroutine destroy_all(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call changes_destroy(changes%component_exchange)
        call changes_destroy(changes%rotation_tuner)
        call changes_destroy(changes%rotated_orientations)
        call changes_destroy(changes%move_tuner)
        call changes_destroy(changes%moved_positions)
    end subroutine destroy_all

    subroutine create_moved_positions(moved_positions, component_positions, &
        periodic_box, input_data, prefix)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: moved_positions
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_moved_positions(moved_positions, component_positions)
        call construct_moved_positions(moved_positions, periodic_box, component_positions, &
            input_data, prefix)
    end subroutine create_moved_positions

    subroutine allocate_moved_positions(moved_positions, component_positions)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: moved_positions
        class(Abstract_Component_Coordinates), intent(in) :: component_positions

        if (component_has_positions(component_positions)) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Changed_Coordinates :: moved_positions)
        end if
    end subroutine allocate_moved_positions

    subroutine construct_moved_positions(moved_positions, periodic_box, component_positions, &
        input_data, prefix)
        class(Abstract_Changed_Coordinates), intent(inout) :: moved_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: move_delta(:)
        type(Concrete_Tuning_Parameters) :: tuning_parameters

        select type (moved_positions)
            type is (Concrete_Moved_Positions)
                data_field = prefix//"Small Move.delta"
                call input_data%get(data_field, move_delta, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Small Move.increase factor"
                call input_data%get(data_field, tuning_parameters%increase_factor, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Small Move.maximum increase factor"
                call input_data%get(data_field, tuning_parameters%increase_factor_max, &
                    data_found)
                call check_data_found(data_field, data_found)
                deallocate(data_field)
                call moved_positions%construct(periodic_box, component_positions, move_delta, &
                    tuning_parameters)
                deallocate(move_delta)
            type is (Null_Changed_Coordinates)
                call moved_positions%construct()
            class default
                call error_exit("construct_moved_positions: type unknown")
        end select
    end subroutine construct_moved_positions

    subroutine destroy_changed_coordinates(changed_coordinates)
        class(Abstract_Changed_Coordinates), allocatable, intent(inout) :: changed_coordinates

        if (allocated(changed_coordinates)) then
            call changed_coordinates%destroy()
            deallocate(changed_coordinates)
        end if
    end subroutine destroy_changed_coordinates

    subroutine changes_create_move_tuner(move_tuner, moved_positions, input_data, prefix)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: move_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: moved_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_move_tuner(move_tuner, moved_positions)
        call construct_change_tuner(move_tuner, moved_positions, input_data, prefix//"Small Move.")
    end subroutine changes_create_move_tuner

    subroutine allocate_move_tuner(move_tuner, moved_positions)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: move_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: moved_positions

        if (component_can_move(moved_positions) .and. num_tuning_steps > 0) then
            allocate(Concrete_Change_Tuner :: move_tuner)
        else
            allocate(Null_Change_Tuner :: move_tuner)
        end if
    end subroutine allocate_move_tuner

    subroutine construct_change_tuner(change_tuner, changed_coordinates, input_data, prefix)
        class(Abstract_Change_Tuner), intent(inout) :: change_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: changed_coordinates
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Change_Tuner_Parameters) :: tuner_parameters

        select type (change_tuner)
            type is (Concrete_Change_Tuner)
                data_field = prefix//"accumulation period"
                call input_data%get(data_field, tuner_parameters%accumulation_period, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"wanted success ratio"
                call input_data%get(data_field, tuner_parameters%wanted_success_ratio, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"tolerance"
                call input_data%get(data_field, tuner_parameters%tolerance, data_found)
                call check_data_found(data_field, data_found)
        end select
        call change_tuner%construct(changed_coordinates, tuner_parameters)
    end subroutine construct_change_tuner

    subroutine destroy_change_tuner(change_tuner)
        class(Abstract_Change_Tuner), allocatable, intent(inout) :: change_tuner

        if (allocated(change_tuner)) then
            call change_tuner%destroy()
            deallocate(change_tuner)
        end if
    end subroutine destroy_change_tuner

    subroutine create_rotated_orientations(rotated_orientations, &
        component_orientations, input_data, prefix)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_rotated_orientations(rotated_orientations, component_orientations)
        call construct_rotated_orientations(rotated_orientations, component_orientations, &
            input_data, prefix)
    end subroutine create_rotated_orientations

    subroutine allocate_rotated_orientations(rotated_orientations, component_orientations)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations

        if (component_has_orientations(component_orientations)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Changed_Coordinates :: rotated_orientations)
        end if
    end subroutine allocate_rotated_orientations

    subroutine construct_rotated_orientations(rotated_orientations, component_orientations, &
        input_data, prefix)
        class(Abstract_Changed_Coordinates), intent(inout) :: rotated_orientations
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: rotation_delta
        type(Concrete_Tuning_Parameters) :: tuning_parameters

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                data_field = prefix//"Small Rotation.delta"
                call input_data%get(data_field, rotation_delta, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Small Rotation.increase factor"
                call input_data%get(data_field, tuning_parameters%increase_factor, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Small Rotation.maximum increase factor"
                call input_data%get(data_field, tuning_parameters%increase_factor_max, &
                    data_found)
                call check_data_found(data_field, data_found)
                deallocate(data_field)
                call rotated_orientations%construct(component_orientations, rotation_delta, &
                    tuning_parameters)
            type is (Null_Changed_Coordinates)
                call rotated_orientations%construct()
            class default
                call error_exit("construct_rotated_orientations: type unknown")
        end select
    end subroutine construct_rotated_orientations

    subroutine changes_create_rotation_tuner(rotation_tuner, rotated_orientations, &
        input_data, prefix)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: rotated_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_rotation_tuner(rotation_tuner, rotated_orientations)
        call construct_change_tuner(rotation_tuner, rotated_orientations, input_data, &
            prefix//"Small Rotation.")
    end subroutine changes_create_rotation_tuner

    subroutine allocate_rotation_tuner(rotation_tuner, rotated_orientations)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: rotated_orientations

        if (component_can_rotate(rotated_orientations) .and. num_tuning_steps > 0) then
            allocate(Concrete_Change_Tuner :: rotation_tuner)
        else
            allocate(Null_Change_Tuner :: rotation_tuner)
        end if
    end subroutine allocate_rotation_tuner

    subroutine create_component_exchange(component_exchange, component)
        class(Abstract_Component_Exchange), allocatable, intent(out) :: component_exchange
        type(Component_Wrapper), intent(in) :: component

        if (component_can_exchange(component%chemical_potential)) then
            allocate(Concrete_Component_Exchange :: component_exchange)
        else
            allocate(Null_Component_Exchange :: component_exchange)
        end if
        call component_exchange%construct(component)
    end subroutine create_component_exchange

    subroutine destroy_component_exchange(component_exchange)
        class(Abstract_Component_Exchange), allocatable, intent(inout) :: component_exchange

        if (allocated(component_exchange)) then
            call component_exchange%destroy()
            deallocate(component_exchange)
        end if
    end subroutine destroy_component_exchange

end module procedures_changes_factory
