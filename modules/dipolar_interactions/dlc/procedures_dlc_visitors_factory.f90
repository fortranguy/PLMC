module procedures_dlc_visitors_factory

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

    subroutine create(visitors, periodic_boxes, reciprocal_lattices, dipolar_interactions_static)
        class(Abstract_DLC_Visitor), allocatable, intent(out) :: visitors(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattices(:)
        type(Dipolar_Interactions_Static_Wrapper), intent(in) ::dipolar_interactions_static(:)

        integer :: i_box

        select type (structures => dipolar_interactions_static(1)%dlc_structures)
            type is (Concrete_DLC_Structures)
                allocate(Concrete_DLC_Visitor :: visitors(size(dipolar_interactions_static)))
            type is (Null_DLC_Structures)
                allocate(Null_DLC_Visitor :: visitors(size(dipolar_interactions_static)))
            class default
                call error_exit("procedures_dlc_visitors_factory: create: structures type unknown.")
        end select

        do i_box = 1, size(dipolar_interactions_static)
            call visitors(i_box)%construct(periodic_boxes(i_box), reciprocal_lattices(i_box), &
                dipolar_interactions_static(i_box)%dlc_weight, dipolar_interactions_static(i_box)%&
                dlc_structures)
        end do
    end subroutine create

    subroutine destroy(visitors)
        class(Abstract_DLC_Visitor), allocatable, intent(inout) :: visitors(:)

        integer :: i_box

        if (allocated(visitors)) then
            do i_box = size(visitors), 1, -1
                call visitors(i_box)%destroy()
            end do
            deallocate(visitors)
        end if
    end subroutine destroy

end module procedures_dlc_visitors_factory
