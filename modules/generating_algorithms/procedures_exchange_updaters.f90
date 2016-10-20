module procedures_exchange_updaters

use types_temporary_particle, only: Concrete_Temporary_Particle
use types_component_wrapper, only: Component_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment
use types_cells_wrapper, only: Cells_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: update_add, update_remove

contains

    subroutine update_add(component, total_moment, cells, dipolar_interactions_static, i_component,&
        particle)
        type(Component_Wrapper), intent(inout) :: component
        class(Abstract_Mixture_Total_Moment), intent(inout) :: total_moment
        type(Cells_Wrapper), intent(inout) :: cells
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component

        call component%num_particles%set(component%num_particles%get() + 1)
        call component%positions%add(particle%position)
        call component%orientations%add(particle%orientation)
        call total_moment%add(i_component, particle%dipole_moment)

        do j_component = 1, size(cells%visitable_cells, 2)
            call cells%visitable_cells(i_component, j_component)%add(particle)
        end do

        call dipolar_interactions_static%reci_structure%update_add(i_component, particle)
        call dipolar_interactions_static%dlc_structures%update_add(i_component, particle)
    end subroutine update_add

    subroutine update_remove(component, total_moment, cells, dipolar_interactions_static, &
        i_component, particle)
        type(Component_Wrapper), intent(inout) :: component
        class(Abstract_Mixture_Total_Moment), intent(inout) :: total_moment
        type(Cells_Wrapper), intent(inout) :: cells
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component

        call dipolar_interactions_static%dlc_structures%update_remove(i_component, particle)
        call dipolar_interactions_static%reci_structure%update_remove(i_component, particle)

        do j_component = size(cells%visitable_cells, 2), 1, -1
            call cells%visitable_cells(i_component, j_component)%remove(particle)
        end do

        call total_moment%remove(i_component, particle%dipole_moment)
        call component%orientations%remove(particle%i)
        call component%positions%remove(particle%i)
        call component%num_particles%set(component%num_particles%get() - 1)
    end subroutine update_remove

end module procedures_exchange_updaters
