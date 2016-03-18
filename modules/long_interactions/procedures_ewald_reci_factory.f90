module procedures_ewald_reci_factory

use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight, Concrete_Ewald_Reci_Weight, &
    Null_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure, &
    Concrete_Ewald_Reci_Structure, Null_Ewald_Reci_Structure
use class_ewald_reci_visitor, only: Abstract_Ewald_Reci_Visitor, &
    Concrete_Ewald_Reci_Visitor, Null_Ewald_Reci_Visitor

implicit none

private
public :: ewald_reci_create, ewald_reci_destroy

interface ewald_reci_create
    module procedure :: create_visitor
    module procedure :: create_structure
    module procedure :: create_weight
end interface ewald_reci_create

interface ewald_reci_destroy
    module procedure :: destroy_weight
    module procedure :: destroy_structure
    module procedure :: destroy_visitor
end interface ewald_reci_destroy

contains

    subroutine create_visitor(visitor, environment, weight, structure)
        class(Abstract_Ewald_Reci_Visitor), allocatable, intent(out) :: visitor
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), intent(in) :: structure

        select type (structure)
            type is (Concrete_Ewald_Reci_Structure)
                allocate(Concrete_Ewald_Reci_Visitor :: visitor)
            type is (Null_Ewald_Reci_Structure)
                allocate(Null_Ewald_Reci_Visitor :: visitor)
            class default
                call error_exit("ewald_reci_create: create_visitor: structure type unknown.")
        end select
        call visitor%construct(environment%periodic_box, environment%reciprocal_lattice, &
            weight, structure)
    end subroutine create_visitor

    subroutine destroy_visitor(visitor)
        class(Abstract_Ewald_Reci_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_visitor

    subroutine create_structure(structure, environment, components, are_dipolar)
        class(Abstract_Ewald_Reci_Structure), allocatable, intent(out) :: structure
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        if (any(are_dipolar)) then
            allocate(Concrete_Ewald_Reci_Structure :: structure)
        else
            allocate(Null_Ewald_Reci_Structure :: structure)
        end if
        call structure%construct(environment%periodic_box, environment%reciprocal_lattice, &
            components, are_dipolar)
    end subroutine create_structure

    subroutine destroy_structure(structure)
        class(Abstract_Ewald_Reci_Structure), allocatable, intent(inout) :: structure

        if (allocated(structure)) then
            call structure%destroy()
            deallocate(structure)
        end if
    end subroutine destroy_structure

    subroutine create_weight(weight, environment, dipoles_exist, alpha)
        class(Abstract_Ewald_Reci_Weight), allocatable, intent(out) :: weight
        type(Environment_Wrapper), intent(in) :: environment
        logical, intent(in) :: dipoles_exist
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha

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
            deallocate(weight)
        end if
    end subroutine destroy_weight

end module procedures_ewald_reci_factory
