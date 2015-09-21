module class_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_positions, only: Abstract_Particles_Positions
use class_moved_positions, only: Abstract_Moved_Positions, &
    Concrete_Moved_Positions
use class_particles_orientations, only:  Abstract_Particles_Orientations
use class_rotated_orientations, only: Abstract_Rotated_Orientations, &
    Concrete_Rotated_Orientations
use class_particles_exchange, only: Particles_Exchange_Facade
use module_adaptation, only: Concrete_Adaptation_Parameters
use types_particles, only: Particles_Wrapper
use types_changes, only: Changes_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Changes_Factory
    contains
        procedure :: create => Abstract_Changes_Factory_create
        procedure, private, nopass ::  allocate_moved_positions => &
            Abstract_Changes_Factory_allocate_moved_positions
        procedure, private, nopass :: set_moved_positions => &
            Abstract_Changes_Factory_set_moved_positions
        procedure, private, nopass :: allocate_rotated_orientations => &
            Abstract_Changes_Factory_allocate_rotated_orientations
        procedure, private, nopass :: set_rotated_orientations => &
            Abstract_Changes_Factory_set_rotated_orientations
        procedure, nopass :: destroy => Abstract_Changes_Factory_destroy
    end type Abstract_Changes_Factory

contains

    subroutine Abstract_Changes_Factory_create(this, changes, input_data, prefix, particles)
        class(Abstract_Changes_Factory), intent(in) :: this
        type(Changes_Wrapper), intent(out) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        type(Particles_Wrapper), intent(in) :: particles

        call this%allocate_moved_positions(changes)
        call this%set_moved_positions(changes, input_data, prefix, particles%positions)
        call this%allocate_rotated_orientations(changes)
        call this%set_rotated_orientations(changes, input_data, prefix, particles%orientations)
        call changes%particles_exchange%construct(particles)
    end subroutine Abstract_Changes_Factory_create

    subroutine Abstract_Changes_Factory_allocate_moved_positions(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        allocate(Concrete_Moved_Positions :: changes%moved_positions)
    end subroutine Abstract_Changes_Factory_allocate_moved_positions

    subroutine Abstract_Changes_Factory_set_moved_positions(changes, input_data, prefix, positions)
        type(Changes_Wrapper), intent(inout) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Positions), intent(in) :: positions

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta(:)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        data_field = prefix//".Small Move.delta"
        call input_data%get(data_field, moved_delta, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Move.increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Move.maximum increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
        call test_data_found(data_field, data_found)
        call changes%moved_positions%construct(positions, moved_delta, adaptation_parameters)
        deallocate(data_field)
    end subroutine Abstract_Changes_Factory_set_moved_positions

    subroutine Abstract_Changes_Factory_allocate_rotated_orientations(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        allocate(Concrete_Rotated_Orientations :: changes%rotated_orientations)
    end subroutine Abstract_Changes_Factory_allocate_rotated_orientations

    subroutine Abstract_Changes_Factory_set_rotated_orientations(changes, input_data, prefix, &
        orientations)
        type(Changes_Wrapper), intent(inout) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Orientations), intent(in) :: orientations

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        data_field = prefix//".Small Rotation.delta"
        call input_data%get(data_field, moved_delta, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Rotation.increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Rotation.maximum increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
        call test_data_found(data_field, data_found)
        call changes%rotated_orientations%construct(orientations, moved_delta, adaptation_parameters)
        deallocate(data_field)
    end subroutine Abstract_Changes_Factory_set_rotated_orientations

    subroutine Abstract_Changes_Factory_destroy(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call changes%particles_exchange%destroy()
        call changes%rotated_orientations%destroy()
        if (allocated(changes%rotated_orientations)) deallocate(changes%rotated_orientations)
        call changes%moved_positions%destroy()
        if (allocated(changes%moved_positions)) deallocate(changes%moved_positions)
    end subroutine Abstract_Changes_Factory_destroy

end module class_changes_factory
