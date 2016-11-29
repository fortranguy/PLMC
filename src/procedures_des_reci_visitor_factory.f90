module procedures_des_reci_visitor_factory

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

    subroutine create(visitor, periodic_box, reciprocal_lattice, dipolar_interactions_static)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Dipolar_Interactions_Static_Wrapper), intent(in) ::dipolar_interactions_static

        select type (structure => dipolar_interactions_static%reci_structure)
            type is (Concrete_DES_Reci_Structure)
                allocate(Concrete_DES_Reci_Visitor :: visitor)
            type is (Null_DES_Reci_Structure)
                allocate(Null_DES_Reci_Visitor :: visitor)
            class default
                call error_exit("procedures_des_reci_visitor_factory: create: structure"//&
                    " type unknown.")
        end select

        call visitor%construct(periodic_box, dipolar_interactions_static%box_size_memento_reci, &
            reciprocal_lattice, dipolar_interactions_static%reci_weight, &
            dipolar_interactions_static%reci_structure)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_DES_Reci_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_des_reci_visitor_factory
