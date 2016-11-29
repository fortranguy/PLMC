module procedures_dlc_visitor_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use classes_dlc_structures, only: Concrete_DLC_Structures, Null_DLC_Structures
use classes_dlc_visitor, only: Abstract_DLC_Visitor, Concrete_DLC_Visitor, Null_DLC_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, periodic_box, reciprocal_lattice, dipolar_interactions_static)
        class(Abstract_DLC_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Dipolar_Interactions_Static_Wrapper), intent(in) ::dipolar_interactions_static

        select type (structures => dipolar_interactions_static%dlc_structures)
            type is (Concrete_DLC_Structures)
                allocate(Concrete_DLC_Visitor :: visitor)
            type is (Null_DLC_Structures)
                allocate(Null_DLC_Visitor :: visitor)
            class default
                call error_exit("procedures_dlc_visitor_factory: create: structures type unknown.")
        end select

        call visitor%construct(periodic_box, reciprocal_lattice, dipolar_interactions_static%&
            dlc_weight, dipolar_interactions_static%dlc_structures)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_DLC_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_dlc_visitor_factory
