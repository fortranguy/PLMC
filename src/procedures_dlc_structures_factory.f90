module procedures_dlc_structures_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use procedures_environment_inquirers, only: periodicity_is_xy
use types_component_wrapper, only: Component_Wrapper
use classes_dlc_structures, only: Abstract_DLC_Structures, Concrete_DLC_Structures, &
    Null_DLC_Structures

implicit none

private
public :: create, destroy

contains

    subroutine create(structures, periodic_box, reciprocal_lattice, components, are_dipolar)
        class(Abstract_DLC_Structures), allocatable, intent(out) :: structures
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (periodicity_is_xy(periodic_box) .and. any(are_dipolar)) then
            allocate(Concrete_DLC_Structures :: structures)
        else
            allocate(Null_DLC_Structures :: structures)
        end if
        call structures%construct(periodic_box, reciprocal_lattice, components, are_dipolar)
    end subroutine create

    subroutine destroy(structures)
        class(Abstract_DLC_Structures), allocatable, intent(inout) :: structures

        if (allocated(structures)) then
            call structures%destroy()
            deallocate(structures)
        end if
    end subroutine destroy

end module procedures_dlc_structures_factory
