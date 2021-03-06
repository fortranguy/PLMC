module procedures_des_real_factory

use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use procedures_des_real_pair_factory, only: des_real_pair_create => create, &
    des_real_pair_destroy => destroy
use classes_des_real_component, only: DES_Real_Component_Wrapper
use procedures_des_real_component_factory, only: des_real_component_create => create, &
    des_real_component_destroy => destroy
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_components
    module procedure :: des_real_component_create
    module procedure :: des_real_pair_create
end interface create

interface destroy
    module procedure :: des_real_pair_destroy
    module procedure :: des_real_component_destroy
    module procedure :: destroy_components
end interface destroy

contains

    subroutine create_components(components, periodic_box, mixture_components, are_dipolar, &
        dipolar_interactions_static)
        type(DES_Real_Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        logical, intent(in) :: are_dipolar(:)
        type(Dipolar_Interactions_Static_Wrapper), intent(in) ::dipolar_interactions_static

        integer :: i_component, j_component
        logical :: interact_ij

        allocate(components(size(mixture_components), size(mixture_components)))

        do j_component = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                interact_ij = are_dipolar(i_component) .and. are_dipolar(j_component)
                call create(components(i_component, j_component)%component, periodic_box, &
                    dipolar_interactions_static%box_size_memento_real, &
                    mixture_components(i_component)%positions, &
                    mixture_components(i_component)%dipole_moments, interact_ij, &
                    dipolar_interactions_static%real_pair)
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(DES_Real_Component_Wrapper), allocatable, intent(inout) :: components(:, :)

        integer :: i_component, j_component

        if (allocated(components)) then
            do j_component = size(components, 2), 1, -1
                do i_component = size(components, 1), 1, -1
                    call destroy(components(i_component, j_component)%component)
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

end module procedures_des_real_factory
