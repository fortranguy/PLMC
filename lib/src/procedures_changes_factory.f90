module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_chemical_potential, only: Abstract_Particles_Chemical_Potential
use class_moved_positions, only: Abstract_Moved_Positions, &
    Concrete_Moved_Positions, Null_Moved_Positions
use class_particles_orientations, only:  Abstract_Particles_Orientations
use class_rotated_orientations, only: Abstract_Rotated_Orientations, &
    Concrete_Rotated_Orientations, Null_Rotated_Orientations
use class_particles_exchange, only: Abstract_Particles_Exchange, &
    Concrete_Particles_Exchange, Null_Particles_Exchange
use module_adaptation, only: Concrete_Adaptation_Parameters
use types_particles, only: Particles_Wrapper
use procedures_types_selectors, only: particles_have_positions, particles_have_orientations, &
    particles_can_exchange
use types_changes, only: Changes_Wrapper

implicit none

private
public :: changes_factory_create, changes_factory_destroy

interface changes_factory_create
    module procedure :: changes_factory_create_all
    module procedure :: allocate_and_construct_moved_positions
    module procedure :: allocate_and_construct_rotated_orientations
    module procedure :: allocate_and_construct_particles_exchange
end interface changes_factory_create

interface changes_factory_destroy
    module procedure :: destroy_and_deallocate_particles_exchange
    module procedure :: destroy_and_deallocate_rotated_orientations
    module procedure :: destroy_and_deallocate_moved_positions
    module procedure :: changes_factory_destroy_all
end interface changes_factory_destroy

contains

    subroutine changes_factory_create_all(changes, input_data, prefix, particles)
        type(Changes_Wrapper), intent(out) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        type(Particles_Wrapper), intent(in) :: particles

        call changes_factory_create(changes%moved_positions, particles%positions, input_data, &
            prefix)
        call changes_factory_create(changes%rotated_orientations, particles%orientations, &
            input_data, prefix)
        call changes_factory_create(changes%particles_exchange, particles)
    end subroutine changes_factory_create_all

    subroutine allocate_and_construct_moved_positions(moved_positions, particles_positions, &
        input_data, prefix)
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_moved_positions(moved_positions, particles_positions)
        call construct_moved_positions(moved_positions, input_data, prefix, particles_positions)
    end subroutine allocate_and_construct_moved_positions

    subroutine allocate_moved_positions(moved_positions, particles_positions)
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions
        class(Abstract_Particles_Positions), intent(in) :: particles_positions

        if (particles_have_positions(particles_positions)) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Moved_Positions :: moved_positions)
        end if
    end subroutine allocate_moved_positions

    subroutine construct_moved_positions(moved_positions, input_data, prefix, particles_positions)
        class(Abstract_Moved_Positions), intent(inout) :: moved_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Positions), intent(in) :: particles_positions

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta(:)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        select type (moved_positions)
            type is (Concrete_Moved_Positions)
                data_field = prefix//"Small Move.delta"
                call input_data%get(data_field, moved_delta, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//"Small Move.increase factor"
                call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//"Small Move.maximum increase factor"
                call input_data%get(data_field, adaptation_parameters%increase_factor_max, &
                    data_found)
                call test_data_found(data_field, data_found)
                deallocate(data_field)
        end select
        call moved_positions%construct(particles_positions, moved_delta, adaptation_parameters)
        if (allocated(moved_delta)) deallocate(moved_delta)
    end subroutine construct_moved_positions

    subroutine allocate_and_construct_rotated_orientations(rotated_orientations, &
        particles_orientations, input_data, prefix)
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_rotated_orientations(rotated_orientations, particles_orientations)
        call construct_rotated_orientations(rotated_orientations, input_data, prefix, &
            particles_orientations)
    end subroutine allocate_and_construct_rotated_orientations

    subroutine allocate_rotated_orientations(rotated_orientations, particles_orientations)
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations

        if (particles_have_orientations(particles_orientations)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Rotated_Orientations :: rotated_orientations)
        end if
    end subroutine allocate_rotated_orientations

    subroutine construct_rotated_orientations(rotated_orientations, input_data, prefix, &
        particles_orientations)
        class(Abstract_Rotated_Orientations), intent(inout) :: rotated_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moved_delta
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                data_field = prefix//"Small Rotation.delta"
                call input_data%get(data_field, moved_delta, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//"Small Rotation.increase factor"
                call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//"Small Rotation.maximum increase factor"
                call input_data%get(data_field, adaptation_parameters%increase_factor_max, &
                    data_found)
                call test_data_found(data_field, data_found)
                deallocate(data_field)
        end select
        call rotated_orientations%construct(particles_orientations, moved_delta, &
            adaptation_parameters)
    end subroutine construct_rotated_orientations

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
        call changes_factory_destroy(changes%rotated_orientations)
        call changes_factory_destroy(changes%moved_positions)
    end subroutine changes_factory_destroy_all

    subroutine destroy_and_deallocate_particles_exchange(particles_exchange)
        class(Abstract_Particles_Exchange), allocatable, intent(inout) :: particles_exchange

        call particles_exchange%destroy()
        if (allocated(particles_exchange)) deallocate(particles_exchange)
    end subroutine destroy_and_deallocate_particles_exchange

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
