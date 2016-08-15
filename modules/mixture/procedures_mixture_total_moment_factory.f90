module procedures_mixture_total_moment_factory

use types_component_wrapper, only: Component_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use procedures_property_inquirers, only: component_is_dipolar

implicit none

private
public :: create, destroy, set_are_dipolar

contains

    subroutine create(total_moment, components)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(out) :: total_moment
        type(Component_Wrapper), intent(in) :: components(:)

        logical :: are_dipolar(size(components))

        call set_are_dipolar(are_dipolar, components)
        if (any(are_dipolar)) then
            allocate(Concrete_Mixture_Total_Moment :: total_moment)
        else
            allocate(Null_Mixture_Total_Moment :: total_moment)
        end if
        call total_moment%construct(components, are_dipolar)
    end subroutine create

    subroutine destroy(total_moment)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(inout) :: total_moment

        if (allocated(total_moment)) then
            call total_moment%destroy()
            deallocate(total_moment)
        end if
    end subroutine destroy

    subroutine set_are_dipolar(are_dipolar, components)
        logical, intent(out) :: are_dipolar(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(are_dipolar)
            are_dipolar(i_component) = component_is_dipolar(components(i_component)%dipole_moments)
        end do
    end subroutine set_are_dipolar

end module procedures_mixture_total_moment_factory
