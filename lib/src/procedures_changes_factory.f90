module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_positions, only: Abstract_Particles_Positions
use class_moved_positions, only: Abstract_Moved_Positions, &
    Concrete_Moved_Positions, Null_Moved_Positions
use class_particles_orientations, only:  Abstract_Particles_Orientations
use class_rotated_orientations, only: Abstract_Rotated_Orientations, &
    Concrete_Rotated_Orientations, Null_Rotated_Orientations
use class_particles_exchange, only: Abstract_Particles_Exchange, &
    Concrete_Particles_Exchange, Null_Particles_Exchange
use module_adaptation, only: Concrete_Adaptation_Parameters
use types_particles, only: Particles_Wrapper
use procedures_particles_factory, only: particles_exist, particles_are_dipolar, &
    particles_can_exchange
use types_changes, only: Changes_Wrapper

implicit none

private
public :: changes_factory_construct, changes_factory_destroy

contains

    subroutine changes_factory_construct(changes, input_data, prefix, particles)
        type(Changes_Wrapper), intent(out) :: changes
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        type(Particles_Wrapper), intent(in) :: particles

        call allocate_moved_positions(changes%moved_positions, input_data, prefix)
        call allocate_rotated_orientations(changes%rotated_orientations, input_data, prefix)
        call allocate_particles_exchange(changes%particles_exchange, input_data, prefix)

        call construct_moved_positions(changes%moved_positions, input_data, prefix, &
            particles%positions)
        call construct_rotated_orientations(changes%rotated_orientations, input_data, prefix, &
            particles%orientations)
        call changes%particles_exchange%construct(particles)
    end subroutine changes_factory_construct

    subroutine allocate_moved_positions(moved_positions, input_data, prefix)
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (particles_exist(input_data, prefix)) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Moved_Positions :: moved_positions)
        end if
    end subroutine allocate_moved_positions

    subroutine allocate_rotated_orientations(rotated_orientations, input_data, prefix)
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (particles_are_dipolar(input_data, prefix)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Rotated_Orientations :: rotated_orientations)
        end if
    end subroutine allocate_rotated_orientations

    subroutine allocate_particles_exchange(particles_exchange, input_data, prefix)
        class(Abstract_Particles_Exchange), allocatable, intent(out) :: particles_exchange
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (particles_can_exchange(input_data, prefix)) then
            allocate(Concrete_Particles_Exchange :: particles_exchange)
        else
            allocate(Null_Particles_Exchange :: particles_exchange)
        end if
    end subroutine allocate_particles_exchange

    subroutine construct_moved_positions(moved_positions, input_data, prefix, positions)
        class(Abstract_Moved_Positions), intent(inout) :: moved_positions
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Positions), intent(in) :: positions

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta(:)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        if (particles_exist(input_data, prefix)) then
            data_field = prefix//".Small Move.delta"
            call input_data%get(data_field, moved_delta, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".Small Move.increase factor"
            call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".Small Move.maximum increase factor"
            call input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
            call test_data_found(data_field, data_found)
            deallocate(data_field)
        end if
        call moved_positions%construct(positions, moved_delta, adaptation_parameters)
        if (allocated(moved_delta)) deallocate(moved_delta)
    end subroutine construct_moved_positions

    subroutine construct_rotated_orientations(rotated_orientations, input_data, prefix, &
        orientations)
        class(Abstract_Rotated_Orientations), intent(inout) :: rotated_orientations
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Orientations), intent(in) :: orientations

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moved_delta
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        if (particles_are_dipolar(input_data, prefix)) then
            data_field = prefix//".Small Rotation.delta"
            call input_data%get(data_field, moved_delta, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".Small Rotation.increase factor"
            call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".Small Rotation.maximum increase factor"
            call input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
            call test_data_found(data_field, data_found)
            deallocate(data_field)
        end if
        call rotated_orientations%construct(orientations, moved_delta, adaptation_parameters)
    end subroutine construct_rotated_orientations

    subroutine changes_factory_destroy(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call changes%particles_exchange%destroy()
        if (allocated(changes%particles_exchange)) deallocate(changes%particles_exchange)
        call changes%rotated_orientations%destroy()
        if (allocated(changes%rotated_orientations)) deallocate(changes%rotated_orientations)
        call changes%moved_positions%destroy()
        if (allocated(changes%moved_positions)) deallocate(changes%moved_positions)
    end subroutine changes_factory_destroy

end module procedures_changes_factory
