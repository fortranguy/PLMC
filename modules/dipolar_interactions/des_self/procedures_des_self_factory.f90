module procedures_des_self_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use types_component_wrapper, only: Component_Wrapper
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use procedures_des_self_component_factory, only: des_self_component_create => create, &
    des_self_component_destroy => destroy
use types_des_self_component_wrapper, only: DES_Self_Component_Wrapper

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

    subroutine create_components(components, periodic_boxes, permittivity, mixture_components, &
        are_dipolar, alpha)
        type(DES_Self_Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) ::periodic_boxes(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Component_Wrapper), intent(in) :: mixture_components(:, :)
        logical, intent(in) :: are_dipolar(:, :)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        integer :: i_box
        integer :: i_component

        allocate(components(size(mixture_components, 1), size(mixture_components, 2)))
        do i_box = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                call create(components(i_component, i_box)%component, periodic_boxes(i_box), &
                    permittivity, mixture_components(i_component, i_box)%dipole_moments, &
                    are_dipolar(i_component, i_box), alpha)
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(DES_Self_Component_Wrapper), allocatable, intent(inout) :: components(:, :)

        integer :: i_box
        integer :: i_component

        if (allocated(components)) then
            do i_box = size(components, 2), 1, -1
                do i_component = size(components, 1), 1, -1
                    call destroy(components(i_component, i_box)%component)
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

end module procedures_des_self_factory
