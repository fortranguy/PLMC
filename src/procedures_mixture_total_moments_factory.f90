module procedures_mixture_total_moments_factory

use types_component_wrapper, only: Component_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use procedures_mixture_inquirers, only: component_is_dipolar

implicit none

private
public :: create, destroy, set_are_dipolar

contains

    !> @todo More restrictive condition for are_dipolar than any?
    subroutine create(total_moments, components)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(out) :: total_moments(:)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box
        logical :: are_dipolar(size(components, 1), size(components, 2))

        call set_are_dipolar(are_dipolar, components)
        if (any(are_dipolar)) then
            allocate(Concrete_Mixture_Total_Moment :: total_moments(size(components, 2)))
        else
            allocate(Null_Mixture_Total_Moment :: total_moments(size(components, 2)))
        end if

        do i_box = 1, size(total_moments)
            call total_moments(i_box)%construct(components(:, i_box), are_dipolar(:, i_box))
        end do
    end subroutine create

    subroutine destroy(total_moments)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(inout) :: total_moments(:)

        integer :: i_box

        if (allocated(total_moments)) then
            do i_box = size(total_moments), 1, -1
                call total_moments(i_box)%destroy()
            end do
            deallocate(total_moments)
        end if
    end subroutine destroy

    subroutine set_are_dipolar(are_dipolar, components)
        logical, intent(inout) :: are_dipolar(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(are_dipolar, 2)
            do i_component = 1, size(are_dipolar, 1)
                are_dipolar(i_component, i_box) = &
                    component_is_dipolar(components(i_component, i_box)%dipole_moments)
            end do
        end do
    end subroutine set_are_dipolar

end module procedures_mixture_total_moments_factory
