module class_changes_factory

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
use module_particles, only: Particles_Wrapper_Parameters, Particles_Wrapper
use types_changes, only: Changes_Wrapper

implicit none

private

    type, public :: Concrete_Changes_Factory
    private
        type(Particles_Wrapper_Parameters) :: parameters
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
    contains
        procedure :: create => Concrete_Changes_Factory_create
        procedure, private ::  allocate_moved_positions => &
            Concrete_Changes_Factory_allocate_moved_positions
        procedure, private :: allocate_rotated_orientations => &
            Concrete_Changes_Factory_allocate_rotated_orientations
        procedure, private :: allocate_particles_exchange => &
            Concrete_Changes_Factory_allocate_particles_exchange
        procedure :: construct => Concrete_Changes_Factory_construct
        procedure, private :: construct_moved_positions => &
            Concrete_Changes_Factory_construct_moved_positions
        procedure, private :: construct_rotated_orientations => &
            Concrete_Changes_Factory_construct_rotated_orientations
        procedure, nopass :: destroy => Concrete_Changes_Factory_destroy
    end type Concrete_Changes_Factory

contains

    subroutine Concrete_Changes_Factory_create(this, changes, parameters, input_data, prefix)
        class(Concrete_Changes_Factory), intent(out) :: this
        type(Changes_Wrapper), intent(out) :: changes
        type(Particles_Wrapper_Parameters) :: parameters
        type(json_file), target, intent(in) :: input_data
        character(len=*), intent(in) :: prefix

        this%parameters = parameters
        this%input_data => input_data
        this%prefix = prefix
        call this%allocate_moved_positions(changes%moved_positions)
        call this%allocate_rotated_orientations(changes%rotated_orientations)
        call this%allocate_particles_exchange(changes%particles_exchange)
    end subroutine Concrete_Changes_Factory_create

    subroutine Concrete_Changes_Factory_allocate_moved_positions(this, moved_positions)
        class(Concrete_Changes_Factory), intent(in) :: this
        class(Abstract_Moved_Positions), allocatable, intent(out) :: moved_positions

        if (this%parameters%exist) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Moved_Positions :: moved_positions)
        end if
    end subroutine Concrete_Changes_Factory_allocate_moved_positions

    subroutine Concrete_Changes_Factory_allocate_rotated_orientations(this, rotated_orientations)
        class(Concrete_Changes_Factory), intent(in) :: this
        class(Abstract_Rotated_Orientations), allocatable, intent(out) :: rotated_orientations

        if (this%parameters%exist .and. this%parameters%are_dipolar) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Rotated_Orientations :: rotated_orientations)
        end if
    end subroutine Concrete_Changes_Factory_allocate_rotated_orientations

    subroutine Concrete_Changes_Factory_allocate_particles_exchange(this, particles_exchange)
        class(Concrete_Changes_Factory), intent(in) :: this
        class(Abstract_Particles_Exchange), allocatable, intent(out) :: particles_exchange

        if (this%parameters%exist .and. this%parameters%can_exchange) then
            allocate(Concrete_Particles_Exchange :: particles_exchange)
        else
            allocate(Null_Particles_Exchange :: particles_exchange)
        end if
    end subroutine Concrete_Changes_Factory_allocate_particles_exchange

    subroutine Concrete_Changes_Factory_construct(this, changes, particles)
        type(Changes_Wrapper), intent(inout) :: changes
        class(Concrete_Changes_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(in) :: particles

        call this%construct_moved_positions(changes, particles%positions)
        call this%construct_rotated_orientations(changes, particles%orientations)
        call changes%particles_exchange%construct(particles)
    end subroutine Concrete_Changes_Factory_construct

    subroutine Concrete_Changes_Factory_construct_moved_positions(this, changes, positions)
        class(Concrete_Changes_Factory), intent(in) :: this
        type(Changes_Wrapper), intent(inout) :: changes
        class(Abstract_Particles_Positions), intent(in) :: positions

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta(:)
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        if (.not. this%parameters%exist) return
        data_field = this%prefix//".Small Move.delta"
        call this%input_data%get(data_field, moved_delta, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Small Move.increase factor"
        call this%input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Small Move.maximum increase factor"
        call this%input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
        call test_data_found(data_field, data_found)
        call changes%moved_positions%construct(positions, moved_delta, adaptation_parameters)
        deallocate(data_field)
    end subroutine Concrete_Changes_Factory_construct_moved_positions

    subroutine Concrete_Changes_Factory_construct_rotated_orientations(this, changes, orientations)
        class(Concrete_Changes_Factory), intent(in) :: this
        type(Changes_Wrapper), intent(inout) :: changes
        class(Abstract_Particles_Orientations), intent(in) :: orientations

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        if (.not. (this%parameters%exist .and. this%parameters%are_dipolar)) return
        data_field = this%prefix//".Small Rotation.delta"
        call this%input_data%get(data_field, moved_delta, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Small Rotation.increase factor"
        call this%input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Small Rotation.maximum increase factor"
        call this%input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
        call test_data_found(data_field, data_found)
        call changes%rotated_orientations%construct(orientations, moved_delta, &
            adaptation_parameters)
        deallocate(data_field)
    end subroutine Concrete_Changes_Factory_construct_rotated_orientations

    subroutine Concrete_Changes_Factory_destroy(this, changes)
        class(Concrete_Changes_Factory), intent(inout) :: this
        type(Changes_Wrapper), intent(inout) :: changes

        call changes%particles_exchange%destroy()
        call changes%rotated_orientations%destroy()
        if (allocated(changes%rotated_orientations)) deallocate(changes%rotated_orientations)
        call changes%moved_positions%destroy()
        if (allocated(changes%moved_positions)) deallocate(changes%moved_positions)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Concrete_Changes_Factory_destroy

end module class_changes_factory
