module procedures_des_reci_structure_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_memento, only: Abstract_Box_Size_Memento
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_component_wrapper, only: Component_Wrapper
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure, Concrete_DES_Reci_Structure, &
    Null_DES_Reci_Structure

implicit none

private
public :: create, destroy

contains

    subroutine create(structure, periodic_box, box_size_memento, reciprocal_lattice, components, &
        are_dipolar)
        class(Abstract_DES_Reci_Structure), allocatable, intent(out) :: structure
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Box_Size_Memento), intent(in) ::box_size_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (any(are_dipolar)) then
            allocate(Concrete_DES_Reci_Structure :: structure)
        else
            allocate(Null_DES_Reci_Structure :: structure)
        end if
        call structure%construct(periodic_box, box_size_memento, reciprocal_lattice, components, &
            are_dipolar)
    end subroutine create

    subroutine destroy(structure)
        class(Abstract_DES_Reci_Structure), allocatable, intent(inout) :: structure

        if (allocated(structure)) then
            call structure%destroy()
            deallocate(structure)
        end if
    end subroutine destroy

end module procedures_des_reci_structure_factory
