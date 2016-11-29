module procedures_transmutation_updaters

use types_particle_wrapper, only: Concrete_Particle
use types_component_wrapper, only: Component_Wrapper
use classes_mixture_total_moment, only: Abstract_Mixture_Total_Moment
use types_cells_wrapper, only: Cells_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private
public :: update

contains

    subroutine update(components, total_moment, cells, dipolar_interactions_static, &
        ij_components, particles)
        type(Component_Wrapper), intent(inout) :: components(:)
        class(Abstract_Mixture_Total_Moment), intent(inout) :: total_moment
        type(Cells_Wrapper), intent(inout) :: cells
        type(Dipolar_Interactions_Static_Wrapper), intent(inout) :: dipolar_interactions_static
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) ::  particles(:)

        integer :: k_component

        do k_component = size(cells%visitable_cells, 2), 1, -1
            call cells%visitable_cells(ij_components(1), k_component)%remove(particles(1))
        end do

        call total_moment%remove(ij_components(1), particles(1)%dipole_moment)
        call components(ij_components(1))%orientations%remove(particles(1)%i)
        call components(ij_components(1))%positions%remove(particles(1)%i)
        call components(ij_components(1))%num_particles%set(components(ij_components(1))%&
            num_particles%get() - 1)

        call components(ij_components(2))%num_particles%set(components(ij_components(2))%&
            num_particles%get() + 1)
        call components(ij_components(2))%positions%add(particles(2)%position)
        call components(ij_components(2))%orientations%add(particles(2)%orientation)
        call total_moment%add(ij_components(2), particles(2)%dipole_moment)

        do k_component = 1, size(cells%visitable_cells, 2)
            call cells%visitable_cells(ij_components(2), k_component)%add(particles(2))
        end do

        call dipolar_interactions_static%reci_structure%&
            update_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
        call dipolar_interactions_static%dlc_structures%&
            update_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
    end subroutine update

end module procedures_transmutation_updaters
