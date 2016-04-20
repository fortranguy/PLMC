module procedures_des_self_component_factory

use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_self_component, only: Abstract_DES_Self_Component, &
    Concrete_DES_Self_Component, Null_DES_Self_Component

implicit none

private
public :: des_self_component_create, des_self_component_destroy

contains

    subroutine des_self_component_create(component, permittivity, dipolar_moments, is_dipolar, alpha)
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
    end subroutine des_self_component_create

    subroutine des_self_component_destroy(component)
        class(Abstract_DES_Self_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine des_self_component_destroy

end module procedures_des_self_component_factory
