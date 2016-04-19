module procedures_change_tuner_factory

use classes_changed_coordinates, only: Abstract_Changed_Coordinates
use types_change_tuner_parameters, only: Concrete_Change_Tuner_Parameters
use classes_change_tuner, only: Abstract_Change_Tuner, Concrete_Change_Tuner, Null_Change_Tuner
use procedures_property_inquirers, only: component_can_move, component_can_rotate
use module_plmc_iterations, only: num_tuning_steps

implicit none

private
public :: change_tuner_create_move, change_tuner_create_rotation, change_tuner_destroy

interface change_tuner_destroy
    module procedure :: destroy_change_tuner
end interface change_tuner_destroy

contains

    subroutine change_tuner_create_move(move_tuner, moved_positions, tuner_parameters)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: move_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: moved_positions
        type(Concrete_Change_Tuner_Parameters), intent(in) :: tuner_parameters

        if (component_can_move(moved_positions) .and. num_tuning_steps > 0) then
            allocate(Concrete_Change_Tuner :: move_tuner)
        else
            allocate(Null_Change_Tuner :: move_tuner)
        end if
        call move_tuner%construct(moved_positions, tuner_parameters)
    end subroutine change_tuner_create_move

    subroutine change_tuner_create_rotation(rotation_tuner, rotated_orientations, tuner_parameters)
        class(Abstract_Change_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Changed_Coordinates), intent(in) :: rotated_orientations
        type(Concrete_Change_Tuner_Parameters), intent(in) :: tuner_parameters

        if (component_can_rotate(rotated_orientations) .and. num_tuning_steps > 0) then
            allocate(Concrete_Change_Tuner :: rotation_tuner)
        else
            allocate(Null_Change_Tuner :: rotation_tuner)
        end if
        call rotation_tuner%construct(rotated_orientations, tuner_parameters)
    end subroutine change_tuner_create_rotation

    subroutine destroy_change_tuner(change_tuner)
        class(Abstract_Change_Tuner), allocatable, intent(inout) :: change_tuner

        if (allocated(change_tuner)) then
            call change_tuner%destroy()
            deallocate(change_tuner)
        end if
    end subroutine destroy_change_tuner

end module procedures_change_tuner_factory
