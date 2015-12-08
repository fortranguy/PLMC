module procedures_ewald_reci_factory

use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight, Concrete_Ewald_Reci_Weight, &
    Null_Ewald_Reci_Weight
use class_ewald_reci_component, only: Abstract_Ewald_Reci_Component, &
    Concrete_Ewald_Reci_Component, Null_Ewald_Reci_Component
use types_long_interactions_wrapper, only: Ewald_Reci_Component_Wrapper

implicit none

private
public :: ewald_reci_create, ewald_reci_destroy

interface ewald_reci_create
    module procedure :: create_components
    module procedure :: create_component
    module procedure :: create_weight
end interface ewald_reci_create

interface ewald_reci_destroy
    module procedure :: destroy_weight
    module procedure :: destroy_component
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
            call ewald_reci_create(reci_components(i_component)%reci_component, environment, &
                components(i_component), are_dipolar(i_component), weight)
        end do
    end subroutine create_components

    subroutine destroy_components(reci_components)
        type(Ewald_Reci_Component_Wrapper), allocatable, intent(inout) :: reci_components(:)

        integer :: i_component

        if (allocated(reci_components)) then
            do i_component = size(reci_components), 1, -1
                call ewald_reci_destroy(reci_components(i_component)%reci_component)
            end do
            deallocate(reci_components)
        end if
    end subroutine destroy_components

    subroutine create_component(reci_component, environment, component, is_dipolar, weight)
        class(Abstract_Ewald_Reci_Component), allocatable, intent(out) :: reci_component
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: component
        logical, intent(in) :: is_dipolar
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight

        if (is_dipolar) then
            allocate(Concrete_Ewald_Reci_Component :: reci_component)
        else
            allocate(Null_Ewald_Reci_Component :: reci_component)
        end if
        call reci_component%construct(environment%periodic_box, environment%reciprocal_lattice, &
            component%positions, component%dipolar_moments, weight)
    end subroutine create_component

    subroutine destroy_component(reci_component)
        class(Abstract_Ewald_Reci_Component), allocatable, intent(inout) :: reci_component

        if (allocated(reci_component)) then
            call reci_component%destroy()
            deallocate(reci_component)
        end if
    end subroutine destroy_component

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
