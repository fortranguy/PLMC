module procedures_des_reci_structure_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use classes_des_reci_structure, only: Abstract_DES_Reci_Structure, Concrete_DES_Reci_Structure, &
    Null_DES_Reci_Structure

implicit none

private
public :: des_reci_structure_create, des_reci_structure_destroy

contains

    subroutine des_reci_structure_create(structure, environment, components, are_dipolar)
        class(Abstract_DES_Reci_Structure), allocatable, intent(out) :: structure
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (any(are_dipolar)) then
            allocate(Concrete_DES_Reci_Structure :: structure)
        else
            allocate(Null_DES_Reci_Structure :: structure)
        end if
        call structure%construct(environment%periodic_box, environment%reciprocal_lattice, &
            components, are_dipolar)
    end subroutine des_reci_structure_create

    subroutine des_reci_structure_destroy(structure)
        class(Abstract_DES_Reci_Structure), allocatable, intent(inout) :: structure

        if (allocated(structure)) then
            call structure%destroy()
            deallocate(structure)
        end if
    end subroutine des_reci_structure_destroy

end module procedures_des_reci_structure_factory
