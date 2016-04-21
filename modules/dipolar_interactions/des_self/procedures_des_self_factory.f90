module procedures_des_self_factory

use classes_permittivity, only: Abstract_Permittivity
use types_component_wrapper, only: Component_Wrapper
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use procedures_des_self_component_factory, only: des_self_component_create => create, &
    des_self_component_destroy => destroy
use types_dipolar_interactions_wrapper, only: DES_Self_Component_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_components
    module procedure :: des_self_component_create
end interface create

interface destroy
    module procedure :: des_self_component_destroy
    module procedure :: destroy_components
end interface destroy

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
            call create(components(i_component)%component, permittivity, &
                mixture_components(i_component)%dipolar_moments, are_dipolar(i_component), alpha)
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(DES_Self_Component_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call destroy(components(i_component)%component)
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

end module procedures_des_self_factory
