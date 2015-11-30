module procedures_ewald_reci_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight, Concrete_Ewald_Reci_Weight, &
    Null_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure, &
    Concrete_Ewald_Reci_Structure, Null_Ewald_Reci_Structure

implicit none

private
public :: ewald_reci_create, ewald_reci_destroy

interface ewald_reci_create
    module procedure :: create_structure
    module procedure :: create_weight
end interface ewald_reci_create

interface ewald_reci_destroy
    module procedure :: destroy_weight
    module procedure :: destroy_structure
end interface ewald_reci_destroy

contains

    subroutine create_structure(structure, environment, component, weight, is_dipolar)
        class(Abstract_Ewald_Reci_Structure), allocatable, intent(out) :: structure
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: component
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        logical, intent(in) :: is_dipolar

        if (is_dipolar) then
            allocate(Concrete_Ewald_Reci_Structure :: structure)
        else
            allocate(Null_Ewald_Reci_Structure :: structure)
        end if
        call structure%construct(environment%periodic_box, environment%reciprocal_lattice, &
            component%positions, component%dipolar_moments, weight)
    end subroutine create_structure

    subroutine destroy_structure(structure)
        class(Abstract_Ewald_Reci_Structure), allocatable, intent(inout) :: structure

        if (allocated(structure)) then
            call structure%destroy()
        end if
    end subroutine destroy_structure

    subroutine create_weight(weight, environment, alpha, dipoles_exist)
        class(Abstract_Ewald_Reci_Weight), allocatable, intent(out) :: weight
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            allocate(Concrete_Ewald_Reci_Weight :: weight)
        else
            allocate(Null_Ewald_Reci_Weight :: weight)
        end if
        call weight%construct(environment%periodic_box, environment%permittivity, environment%&
            reciprocal_lattice, alpha)
    end subroutine create_weight

    subroutine destroy_weight(weight)
        class(Abstract_Ewald_Reci_Weight), allocatable, intent(inout) :: weight

        if (allocated(weight)) then
            call weight%destroy()
        end if
    end subroutine destroy_weight

end module procedures_ewald_reci_factory
