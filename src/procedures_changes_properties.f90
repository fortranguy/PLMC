module procedures_changes_properties

use procedures_mixture_inquirers, only: component_can_translate, component_can_rotate
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private
public :: set_can_translate, set_can_rotate

contains

    subroutine set_can_translate(can_translate, components)
        logical, intent(inout) :: can_translate(:, :)
        type(Changes_Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_translate, 2)
            do i_component = 1, size(can_translate, 1)
                can_translate(i_component, i_box) = &
                    component_can_translate(components(i_component, i_box)%translated_positions)
            end do
        end do
    end subroutine set_can_translate

    subroutine set_can_rotate(can_rotate, components)
        logical, intent(inout) :: can_rotate(:, :)
        type(Changes_Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_rotate, 2)
            do i_component = 1, size(can_rotate, 1)
                can_rotate(i_component, i_box) = &
                    component_can_rotate(components(i_component, i_box)%rotated_orientations)
            end do
        end do
    end subroutine set_can_rotate

end module procedures_changes_properties
