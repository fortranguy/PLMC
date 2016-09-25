module procedures_des_self_component_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_self_component, only: Abstract_DES_Self_Component, &
    Concrete_DES_Self_Component, Null_DES_Self_Component

implicit none

private
public :: create, destroy

contains

    subroutine create(component, periodic_box, permittivity, dipole_moments, is_dipolar, alpha)
        class(Abstract_DES_Self_Component), allocatable, intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) ::periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments
        logical, intent(in) :: is_dipolar
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        if (is_dipolar) then
            allocate(Concrete_DES_Self_Component :: component)
        else
            allocate(Null_DES_Self_Component :: component)
        end if
        call component%construct(periodic_box, permittivity, dipole_moments, alpha)
    end subroutine create

    subroutine destroy(component)
        class(Abstract_DES_Self_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy

end module procedures_des_self_component_factory
