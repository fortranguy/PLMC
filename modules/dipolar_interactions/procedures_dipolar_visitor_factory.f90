module procedures_dipolar_visitor_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use procedures_environment_inquirers, only: periodicity_is_xyz
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use classes_dipolar_visitor, only: Abstract_Dipolar_Visitor, Scalable_Dipolar_Visitor, &
    Unscalable_Dipolar_Visitor, Null_Dipolar_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, periodic_box, box_size_can_change, components, dipolar_interactions)
        class(Abstract_Dipolar_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: box_size_can_change
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        logical :: are_dipolar(size(components))

        call set_are_dipolar(are_dipolar, components)
        if (box_size_can_change .and. any(are_dipolar)) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(Scalable_Dipolar_Visitor :: visitor)
            else
                allocate(Unscalable_Dipolar_Visitor :: visitor)
            end if
        else
            allocate(Null_Dipolar_Visitor :: visitor)
        end if

        select type (visitor)
            type is (Scalable_Dipolar_Visitor)
            type is (Unscalable_Dipolar_Visitor)
                call visitor%construct(components, dipolar_interactions)
            type is (Null_Dipolar_Visitor)
            class default
                call error_exit("procedures_dipolar_visitor_factory: create: visitor type is "//&
                    "unknown.")
        end select
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_Dipolar_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_dipolar_visitor_factory
