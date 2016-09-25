module procedures_des_reci_structure_factory

use classes_box_size_memento, only: Abstract_Box_Size_Memento
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure, Concrete_DES_Reci_Structure, &
    Null_DES_Reci_Structure

implicit none

private
public :: create, destroy

contains

    subroutine create(structure, environment, box_size_memento, components, are_dipolar)
        class(Abstract_DES_Reci_Structure), allocatable, intent(out) :: structure
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Box_Size_Memento), intent(in) ::box_size_memento
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (any(are_dipolar)) then
            allocate(Concrete_DES_Reci_Structure :: structure)
        else
            allocate(Null_DES_Reci_Structure :: structure)
        end if
        call structure%construct(environment%periodic_box, box_size_memento, environment%&
            reciprocal_lattice, components, are_dipolar)
    end subroutine create

    subroutine destroy(structure)
        class(Abstract_DES_Reci_Structure), allocatable, intent(inout) :: structure

        if (allocated(structure)) then
            call structure%destroy()
            deallocate(structure)
        end if
    end subroutine destroy

end module procedures_des_reci_structure_factory
