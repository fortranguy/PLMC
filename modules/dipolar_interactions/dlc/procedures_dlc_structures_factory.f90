module procedures_dlc_structures_factory

use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: periodicity_is_xy
use types_component_wrapper, only: Component_Wrapper
use classes_dlc_structures, only: Abstract_DLC_Structures, Concrete_DLC_Structures, &
    Null_DLC_Structures

implicit none

private
public :: create, destroy

contains

    subroutine create(structures, environment, components, are_dipolar)
        class(Abstract_DLC_Structures), allocatable, intent(out) :: structures
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (periodicity_is_xy(environment%periodic_box) .and. any(are_dipolar)) then
            allocate(Concrete_DLC_Structures :: structures)
        else
            allocate(Null_DLC_Structures :: structures)
        end if
        call structures%construct(environment%periodic_box, environment%reciprocal_lattice, &
            components, are_dipolar)
    end subroutine create

    subroutine destroy(structures)
        class(Abstract_DLC_Structures), allocatable, intent(inout) :: structures

        if (allocated(structures)) then
            call structures%destroy()
            deallocate(structures)
        end if
    end subroutine destroy

end module procedures_dlc_structures_factory
