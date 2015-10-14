module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_plmc_iterations, only: num_tuning_steps
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_component_positions, only: Abstract_Component_Positions
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential
use class_changed_coordinates, only: Abstract_Changed_Coordinates
use class_moved_positions, only: Abstract_Moved_Positions, &
    Concrete_Moved_Positions, Null_Moved_Positions
use class_component_orientations, only:  Abstract_Component_Orientations
use class_rotated_orientations, only: Abstract_Rotated_Orientations, &
    Concrete_Rotated_Orientations, Null_Rotated_Orientations
use class_change_tuner, only: Concrete_Change_Tuner_Parameters, &
    Abstract_Change_Tuner, Concrete_Change_Tuner, Null_Change_Tuner
use class_particles_exchange, only: Abstract_Particles_Exchange, &
    Concrete_Particles_Exchange, Null_Particles_Exchange
use module_change_tuning, only: Concrete_Tuning_Parameters
use types_particles_wrapper, only: Particles_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_property_inquirers, only: particles_have_positions, particles_have_orientations, &
    particles_can_move, particles_can_rotate, particles_can_exchange

implicit none

private
public :: changes_factory_create, changes_factory_destroy

interface changes_factory_create
    module procedure :: changes_factory_create_all
    module procedure :: allocate_and_construct_moved_positions
    module procedure :: allocate_and_construct_move_tuner
    module procedure :: allocate_and_construct_rotated_orientations
    module procedure :: allocate_and_construct_rotation_tuner
    module procedure :: allocate_and_construct_particles_exchange
end interface changes_factory_create

interface changes_factory_destroy
    module procedure :: destroy_and_deallocate_particles_exchange
    module procedure :: destroy_and_deallocate_change_tuner
    module procedure :: destroy_and_deallocate_rotated_orientations
    module procedure :: destroy_and_deallocate_moved_positions
    module procedure :: changes_factory_destroy_all
end interface changes_factory_destroy

contains

    subroutine changes_factory_create_all(changes, periodic_box, particles, input_data, prefix)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Particles_Wrapper), intent(in) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call changes_factory_create(changes%moved_positions, particles%positions, periodic_box, &
            input_data, prefix)
        call changes_factory_create(changes%move_tuner, changes%moved_positions, input_data, prefix)
        call changes_factory_create(changes%rotated_orientations, particles%orientations, &
            input_data, prefix)
        call changes_factory_create(changes%rotation_tuner, changes%rotated_orientations, &
            input_data, prefix)
        call changes_factory_create(changes%particles_exchange, particles)
    end subroutine changes_factory_create_all

    subroutine allocate_and_construct_moved_positions(moved_positions, particles_positions, &
        periodic_box, input_data, prefix)
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions
        class(Abstract_Component_Positions), intent(in) :: particles_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_moved_positions(moved_positions, particles_positions)
        call construct_moved_positions(moved_positions, periodic_box, particles_positions, &
            input_data, prefix)
    end subroutine allocate_and_construct_moved_positions

    subroutine allocate_moved_positions(moved_positions, particles_positions)
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions
        class(Abstract_Component_Positions), intent(in) :: particles_positions

        if (particles_have_positions(particles_positions)) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Moved_Positions :: moved_positions)
        end if
    end subroutine allocate_moved_positions

    subroutine construct_moved_positions(moved_positions, periodic_box, particles_positions, &
        input_data, prefix)
        class(Abstract_Moved_Positions), intent(inout) :: moved_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Positions), intent(in) :: particles_positions
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
        end select
        call moved_positions%construct(periodic_box, particles_positions, move_delta, &
            tuning_parameters)
        if (allocated(move_delta)) deallocate(move_delta)
    end subroutine construct_moved_positions

    subroutine allocate_and_construct_move_tuner(move_tuner, moved_positions, input_data, prefix)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: move_tuner
        class(Abstract_Moved_Positions), intent(in) :: moved_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_move_tuner(move_tuner, moved_positions)
        call construct_change_tuner(move_tuner, moved_positions, input_data, prefix//"Small Move.")
    end subroutine allocate_and_construct_move_tuner

    subroutine allocate_move_tuner(move_tuner, moved_positions)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: move_tuner
        class(Abstract_Moved_Positions), intent(in) :: moved_positions

        if (particles_can_move(moved_positions) .and. num_tuning_steps > 0) then
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

    subroutine allocate_and_construct_rotated_orientations(rotated_orientations, &
        particles_orientations, input_data, prefix)
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Component_Orientations), intent(in) :: particles_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_rotated_orientations(rotated_orientations, particles_orientations)
        call construct_rotated_orientations(rotated_orientations, particles_orientations, &
            input_data, prefix)
    end subroutine allocate_and_construct_rotated_orientations

    subroutine allocate_rotated_orientations(rotated_orientations, particles_orientations)
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Component_Orientations), intent(in) :: particles_orientations

        if (particles_have_orientations(particles_orientations)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Rotated_Orientations :: rotated_orientations)
        end if
    end subroutine allocate_rotated_orientations

    subroutine construct_rotated_orientations(rotated_orientations, particles_orientations, &
        input_data, prefix)
        class(Abstract_Rotated_Orientations), intent(inout) :: rotated_orientations
        class(Abstract_Component_Orientations), intent(in) :: particles_orientations
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
        end select
        call rotated_orientations%construct(particles_orientations, rotation_delta, &
            tuning_parameters)
    end subroutine construct_rotated_orientations

    subroutine allocate_and_construct_rotation_tuner(rotation_tuner, rotated_orientations, &
        input_data, prefix)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_rotation_tuner(rotation_tuner, rotated_orientations)
        call construct_change_tuner(rotation_tuner, rotated_orientations, input_data, &
            prefix//"Small Rotation.")
    end subroutine allocate_and_construct_rotation_tuner

    subroutine allocate_rotation_tuner(rotation_tuner, rotated_orientations)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations

        if (particles_can_rotate(rotated_orientations) .and. num_tuning_steps > 0) then
            allocate(Concrete_Change_Tuner :: rotation_tuner)
        else
            allocate(Null_Change_Tuner :: rotation_tuner)
        end if
    end subroutine allocate_rotation_tuner

    subroutine allocate_and_construct_particles_exchange(particles_exchange, particles)
        class(Abstract_Particles_Exchange), allocatable, intent(out) :: particles_exchange
        type(Particles_Wrapper), intent(in) :: particles

        if (particles_can_exchange(particles%chemical_potential)) then
            allocate(Concrete_Particles_Exchange :: particles_exchange)
        else
            allocate(Null_Particles_Exchange :: particles_exchange)
        end if
        call particles_exchange%construct(particles)
    end subroutine allocate_and_construct_particles_exchange

    subroutine changes_factory_destroy_all(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call changes_factory_destroy(changes%particles_exchange)
        call changes_factory_destroy(changes%rotation_tuner)
        call changes_factory_destroy(changes%rotated_orientations)
        call changes_factory_destroy(changes%move_tuner)
        call changes_factory_destroy(changes%moved_positions)
    end subroutine changes_factory_destroy_all

    subroutine destroy_and_deallocate_particles_exchange(particles_exchange)
        class(Abstract_Particles_Exchange), allocatable, intent(inout) :: particles_exchange

        call particles_exchange%destroy()
        if (allocated(particles_exchange)) deallocate(particles_exchange)
    end subroutine destroy_and_deallocate_particles_exchange

    subroutine destroy_and_deallocate_change_tuner(change_tuner)
        class(Abstract_Change_Tuner), allocatable, intent(inout) :: change_tuner

        call change_tuner%destroy()
        if (allocated(change_tuner)) deallocate(change_tuner)
    end subroutine destroy_and_deallocate_change_tuner

    subroutine destroy_and_deallocate_rotated_orientations(rotated_orientations)
        class(Abstract_Rotated_Orientations), allocatable, intent(inout) :: rotated_orientations

        call rotated_orientations%destroy()
        if (allocated(rotated_orientations)) deallocate(rotated_orientations)
    end subroutine destroy_and_deallocate_rotated_orientations

    subroutine destroy_and_deallocate_moved_positions(moved_positions)
        class(Abstract_Moved_Positions), allocatable, intent(inout) :: moved_positions

        call moved_positions%destroy()
        if (allocated(moved_positions)) deallocate(moved_positions)
    end subroutine destroy_and_deallocate_moved_positions

end module procedures_changes_factory
