module procedures_move_tuner_factory

use classes_changed_box_size, only: Abstract_Changed_Box_Size
use procedures_environment_inquirers, only: box_size_can_change
use procedures_mixture_inquirers, only: component_can_translate, component_can_rotate
use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use classes_move_tuner, only: Abstract_Move_Tuner, Concrete_Move_Tuner, Null_Move_Tuner

implicit none

private
public :: create_boxes_size_change, create_translation, create_rotation, destroy

interface destroy
    module procedure :: destroy_line
    module procedure :: destroy_element
end interface destroy

contains

    subroutine create_boxes_size_change(boxes_size_change_tuner, changed_boxes_size, &
        tuner_parameters, num_tuning_steps)
        class(Abstract_Move_Tuner), allocatable, intent(out) :: boxes_size_change_tuner(:)
        class(Abstract_Changed_Box_Size), intent(in) :: changed_boxes_size(:)
        type(Concrete_Move_Tuner_Parameters), intent(in) :: tuner_parameters
        integer, intent(in) :: num_tuning_steps

        integer :: i_box

        if (all(box_size_can_change(changed_boxes_size)) .and. num_tuning_steps > 0) then
            allocate(Concrete_Move_Tuner :: boxes_size_change_tuner(size(changed_boxes_size)))
        else
            allocate(Null_Move_Tuner :: boxes_size_change_tuner(size(changed_boxes_size)))
        end if

        do i_box = 1, size(boxes_size_change_tuner)
            call boxes_size_change_tuner(i_box)%construct(changed_boxes_size(i_box), &
                tuner_parameters, num_tuning_steps)
        end do
    end subroutine create_boxes_size_change

    subroutine create_box_size_change(box_size_change_tuner, changed_box_size, tuner_parameters, &
        num_tuning_steps)
        class(Abstract_Move_Tuner), allocatable, intent(out) :: box_size_change_tuner
        class(Abstract_Changed_Box_Size), intent(in) :: changed_box_size
        type(Concrete_Move_Tuner_Parameters), intent(in) :: tuner_parameters
        integer, intent(in) :: num_tuning_steps

        if (box_size_can_change(changed_box_size) .and. num_tuning_steps > 0) then
            allocate(Concrete_Move_Tuner :: box_size_change_tuner)
        else
            allocate(Null_Move_Tuner :: box_size_change_tuner)
        end if
        call box_size_change_tuner%construct(changed_box_size, tuner_parameters, num_tuning_steps)
    end subroutine create_box_size_change

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

    subroutine destroy_line(move_tuners)
        class(Abstract_Move_Tuner), allocatable, intent(inout) :: move_tuners(:)

        integer :: i_element

        if (allocated(move_tuners)) then
            do i_element = size(move_tuners), 1, -1
                call move_tuners(i_element)%destroy()
            end do
            deallocate(move_tuners)
        end if
    end subroutine destroy_line

    subroutine destroy_element(move_tuner)
        class(Abstract_Move_Tuner), allocatable, intent(inout) :: move_tuner

        if (allocated(move_tuner)) then
            call move_tuner%destroy()
            deallocate(move_tuner)
        end if
    end subroutine destroy_element

end module procedures_move_tuner_factory
