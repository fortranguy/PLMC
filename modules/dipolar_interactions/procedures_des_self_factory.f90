module procedures_des_self_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use types_component_wrapper, only: Component_Wrapper
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use classes_des_self_component, only: Abstract_DES_Self_Component, &
    Concrete_DES_Self_Component, Null_DES_Self_Component
use classes_des_self_component, only: DES_Self_Component_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface destroy

contains

    subroutine create_line(components, periodic_box, permittivity, mixture_components, are_dipolar,&
        alpha)
        type(DES_Self_Component_Wrapper), allocatable, intent(out) :: components(:)
        class(Abstract_Periodic_Box), intent(in) ::periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        logical, intent(in) :: are_dipolar(:)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        integer :: i_component

        allocate(components(size(mixture_components)))
        do i_component = 1, size(components)
            call create(components(i_component)%component, periodic_box, &
                permittivity, mixture_components(i_component)%dipole_moments, &
                are_dipolar(i_component), alpha)
        end do
    end subroutine create_line

    subroutine destroy_line(components)
        type(DES_Self_Component_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call destroy(components(i_component)%component)
            end do
            deallocate(components)
        end if
    end subroutine destroy_line

    subroutine create_element(component, periodic_box, permittivity, dipole_moments, is_dipolar, &
        alpha)
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
    end subroutine create_element

    subroutine destroy_element(component)
        class(Abstract_DES_Self_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy_element

end module procedures_des_self_factory
