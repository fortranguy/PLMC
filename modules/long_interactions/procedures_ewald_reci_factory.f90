module procedures_ewald_reci_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight, Concrete_Ewald_Reci_Weight, &
    Null_Ewald_Reci_Weight
use class_ewald_reci_structure_, only: Abstract_Ewald_Reci_Structure_, &
    Concrete_Ewald_Reci_Structure_, Null_Ewald_Reci_Structure_
use class_ewald_reci_delta_visitor, only: Abstract_Ewald_Reci_Delta_Visitor, &
    Concrete_Ewald_Reci_Delta_Visitor, Null_Ewald_Reci_Delta_Visitor
use types_long_interactions_wrapper, only: Ewald_Reci_Component_Wrapper

implicit none

private
public :: ewald_reci_create, ewald_reci_destroy

interface ewald_reci_create
    module procedure :: create_components
    module procedure :: create_delta_visitor
    module procedure :: create_structure
    module procedure :: create_weight
end interface ewald_reci_create

interface ewald_reci_destroy
    module procedure :: destroy_weight
    module procedure :: destroy_structure
    module procedure :: destroy_delta_visitor
    module procedure :: destroy_components
end interface ewald_reci_destroy

contains

    subroutine create_components(reci_components, environment, components, are_dipolar, weight)
        type(Ewald_Reci_Component_Wrapper), allocatable, intent(out) :: reci_components(:)
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight

        integer :: i_component
        allocate(reci_components(size(components)))
        do i_component = 1, size(reci_components)
            call ewald_reci_create(reci_components(i_component)%reci_structure, environment, &
                components(i_component), are_dipolar(i_component))
            call ewald_reci_create(reci_components(i_component)%reci_delta_visitor, &
                environment, are_dipolar(i_component), weight, reci_components(i_component)%&
                reci_structure)
        end do
    end subroutine create_components

    subroutine destroy_components(reci_components)
        type(Ewald_Reci_Component_Wrapper), allocatable, intent(inout) :: reci_components(:)

        integer :: i_component

        if (allocated(reci_components)) then
            do i_component = size(reci_components), 1, -1
                call ewald_reci_destroy(reci_components(i_component)%reci_delta_visitor)
                call ewald_reci_destroy(reci_components(i_component)%reci_structure)
            end do
            deallocate(reci_components)
        end if
    end subroutine destroy_components

    subroutine create_delta_visitor(reci_delta_visitor, environment, is_dipolar, weight, structure)
        class(Abstract_Ewald_Reci_Delta_Visitor), allocatable, intent(out) :: reci_delta_visitor
        type(Environment_Wrapper), intent(in) :: environment
        logical, intent(in) :: is_dipolar
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure_), intent(in) :: structure

        if (is_dipolar) then
            allocate(Concrete_Ewald_Reci_Delta_Visitor :: reci_delta_visitor)
        else
            allocate(Null_Ewald_Reci_Delta_Visitor :: reci_delta_visitor)
        end if
        call reci_delta_visitor%construct(environment%periodic_box, environment%&
            reciprocal_lattice, weight, structure)
    end subroutine create_delta_visitor

    subroutine destroy_delta_visitor(reci_delta_visitor)
        class(Abstract_Ewald_Reci_Delta_Visitor), allocatable, intent(inout) :: reci_delta_visitor

        if (allocated(reci_delta_visitor)) then
            call reci_delta_visitor%destroy()
            deallocate(reci_delta_visitor)
        end if
    end subroutine destroy_delta_visitor

    subroutine create_structure(reci_structure, environment, component, is_dipolar)
        class(Abstract_Ewald_Reci_Structure_), allocatable, intent(out) :: reci_structure
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: component
        logical, intent(in) :: is_dipolar

        if (is_dipolar) then
            allocate(Concrete_Ewald_Reci_Structure_ :: reci_structure)
        else
            allocate(Null_Ewald_Reci_Structure_ :: reci_structure)
        end if
        call reci_structure%construct(environment%periodic_box, environment%reciprocal_lattice, &
            component%positions, component%dipolar_moments)
    end subroutine create_structure

    subroutine destroy_structure(reci_structure)
        class(Abstract_Ewald_Reci_Structure_), allocatable, intent(inout) :: reci_structure

        if (allocated(reci_structure)) then
            call reci_structure%destroy()
            deallocate(reci_structure)
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
