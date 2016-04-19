module procedures_des_self_factory

use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_self_component, only: Abstract_DES_Self_Component, &
    Concrete_DES_Self_Component, Null_DES_Self_Component
use types_dipolar_interactions_wrapper, only: DES_Self_Component_Wrapper

implicit none

private
public :: des_self_create, des_self_destroy

interface des_self_create
    module procedure :: create_components
    module procedure :: create_component
end interface des_self_create

interface des_self_destroy
    module procedure :: destroy_component
    module procedure :: destroy_components
end interface des_self_destroy

contains

    subroutine create_components(components, permittivity, mixture_components, are_dipolar, alpha)
        type(DES_Self_Component_Wrapper), allocatable, intent(out) :: components(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        logical, intent(in) :: are_dipolar(:)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        integer :: i_component

        allocate(components(size(mixture_components)))
        do i_component = 1, size(components)
            call des_self_create(components(i_component)%component, permittivity, &
                mixture_components(i_component)%dipolar_moments, are_dipolar(i_component), alpha)
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(DES_Self_Component_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call des_self_destroy(components(i_component)%component)
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine create_component(component, permittivity, dipolar_moments, is_dipolar, alpha)
        class(Abstract_DES_Self_Component), allocatable, intent(out) :: component
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        logical, intent(in) :: is_dipolar
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        if (is_dipolar) then
            allocate(Concrete_DES_Self_Component :: component)
        else
            allocate(Null_DES_Self_Component :: component)
        end if
        call component%construct(permittivity, dipolar_moments, alpha)
    end subroutine create_component

    subroutine destroy_component(component)
        class(Abstract_DES_Self_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy_component

end module procedures_des_self_factory
