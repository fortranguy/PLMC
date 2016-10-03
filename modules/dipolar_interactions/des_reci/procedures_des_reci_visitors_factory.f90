module procedures_des_reci_visitors_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_des_reci_structure, only: Concrete_DES_Reci_Structure, Null_DES_Reci_Structure
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use classes_des_reci_visitor, only: Abstract_DES_Reci_Visitor, &
    Concrete_DES_Reci_Visitor, Null_DES_Reci_Visitor

implicit none

private
public :: create, destroy

contains

    !> @todo
    !> Better select type
    subroutine create(visitors, periodic_boxes, reciprocal_lattices, dipolar_interactions_static)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(out) :: visitors(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattices(:)
        type(Dipolar_Interactions_Static_Wrapper), intent(in) ::dipolar_interactions_static(:)

        integer :: i_box

        select type (structure => dipolar_interactions_static(1)%reci_structure)
            type is (Concrete_DES_Reci_Structure)
                allocate(Concrete_DES_Reci_Visitor :: visitors(size(dipolar_interactions_static)))
            type is (Null_DES_Reci_Structure)
                allocate(Null_DES_Reci_Visitor :: visitors(size(dipolar_interactions_static)))
            class default
                call error_exit("procedures_des_reci_visitors_factory: create: structure"//&
                    " type unknown.")
        end select

        do i_box = 1, size(visitors)
            call visitors(i_box)%construct(periodic_boxes(i_box), &
                dipolar_interactions_static(i_box)%box_size_memento_reci, &
                reciprocal_lattices(i_box), dipolar_interactions_static(i_box)%reci_weight, &
                dipolar_interactions_static(i_box)%reci_structure)
        end do
    end subroutine create

    subroutine destroy(visitors)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(inout) :: visitors(:)

        integer :: i_box

        if (allocated(visitors)) then
            do i_box = size(visitors), 1, -1
                call visitors(i_box)%destroy()
            end do
            deallocate(visitors)
        end if
    end subroutine destroy

end module procedures_des_reci_visitors_factory
