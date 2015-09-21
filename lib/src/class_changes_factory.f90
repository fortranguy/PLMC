module class_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_moved_positions, only: Abstract_Moved_Positions, &
    Concrete_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_particles_exchange, only: Particles_Exchange_Facade
use module_adaptation, only: Concrete_Adaptation_Parameters
use types_changes, only: Changes_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Changes_Factory
    contains
        procedure :: create => Abstract_Changes_Factory_create
        procedure, private, nopass :: set_moved_positions => &
            Abstract_Changes_Factory_set_moved_positions
        procedure, nopass :: destroy => Abstract_Changes_Factory_destroy
    end type Abstract_Changes_Factory

contains

    subroutine Abstract_Changes_Factory_create(this, changes, input_data, prefix)
        class(Abstract_Changes_Factory), intent(in) :: this
        type(Changes_Wrapper), intent(out) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        !call changes%moved_positions%construct(positions, moved_delta, adaptation_parameters)

    end subroutine Abstract_Changes_Factory_create

    subroutine Abstract_Changes_Factory_set_moved_positions(changes, input_data, prefix)
        type(Changes_Wrapper), intent(inout) :: changes
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: moved_delta
        type(Concrete_Adaptation_Parameters) :: adaptation_parameters

        allocate(Concrete_Moved_Positions :: changes%moved_positions)
        data_field = prefix//".Small Move.delta"
        call input_data%get(data_field, moved_delta, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Move.increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor, data_found)
        call test_data_found(data_field, data_found)
        data_field = prefix//".Small Move.maximum increase factor"
        call input_data%get(data_field, adaptation_parameters%increase_factor_max, data_found)
        call test_data_found(data_field, data_found)
    end subroutine Abstract_Changes_Factory_set_moved_positions

    subroutine Abstract_Changes_Factory_destroy(changes)
        type(Changes_Wrapper), intent(inout) :: changes

    end subroutine Abstract_Changes_Factory_destroy

end module class_changes_factory
