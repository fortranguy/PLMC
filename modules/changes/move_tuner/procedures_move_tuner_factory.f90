module procedures_move_tuner_factory

use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use classes_move_tuner, only: Abstract_Move_Tuner, Concrete_Move_Tuner, Null_Move_Tuner
use procedures_property_inquirers, only: component_can_translate, component_can_rotate

implicit none

private
public :: create_translation, create_rotation, destroy

contains

    subroutine create_translation(translation_tuner, translated_positions, tuner_parameters, &
        num_tuning_steps)
        class(Abstract_Move_Tuner), allocatable, intent(out) :: translation_tuner
        class(Abstract_Moved_Coordinates), intent(in) :: translated_positions
        type(Concrete_Move_Tuner_Parameters), intent(in) :: tuner_parameters
        integer, intent(in) :: num_tuning_steps

        if (component_can_translate(translated_positions) .and. num_tuning_steps > 0) then
            allocate(Concrete_Move_Tuner :: translation_tuner)
        else
            allocate(Null_Move_Tuner :: translation_tuner)
        end if
        call translation_tuner%construct(translated_positions, tuner_parameters, num_tuning_steps)
    end subroutine create_translation

    subroutine create_rotation(rotation_tuner, rotated_orientations, tuner_parameters, &
        num_tuning_steps)
        class(Abstract_Move_Tuner), allocatable, intent(out) :: rotation_tuner
        class(Abstract_Moved_Coordinates), intent(in) :: rotated_orientations
        type(Concrete_Move_Tuner_Parameters), intent(in) :: tuner_parameters
        integer, intent(in) :: num_tuning_steps

        if (component_can_rotate(rotated_orientations) .and. num_tuning_steps > 0) then
            allocate(Concrete_Move_Tuner :: rotation_tuner)
        else
            allocate(Null_Move_Tuner :: rotation_tuner)
        end if
        call rotation_tuner%construct(rotated_orientations, tuner_parameters, num_tuning_steps)
    end subroutine create_rotation

    subroutine destroy(move_tuner)
        class(Abstract_Move_Tuner), allocatable, intent(inout) :: move_tuner

        if (allocated(move_tuner)) then
            call move_tuner%destroy()
            deallocate(move_tuner)
        end if
    end subroutine destroy

end module procedures_move_tuner_factory
