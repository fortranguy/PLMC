module procedures_swap_visitors

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_visit_condition, only: visit_different
use classes_external_field, only: Abstract_External_Field
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_pair_potential, only: Pair_Potential_Wrapper
use types_cells_wrapper, only: Cells_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use procedures_dipoles_field_interaction, only:  dipoles_field_visit_add => visit_add, &
    dipoles_field_visit_remove => visit_remove
use types_particle_wrapper, only: Concrete_Particle

implicit none

private
public :: transmutation_visit_short, transmutation_visit_dipolar, i_exclude_particle

interface transmutation_visit_short
    module procedure :: visit_walls
    module procedure :: visit_short
end interface transmutation_visit_short

interface transmutation_visit_dipolar
    module procedure :: visit_field
    module procedure :: visit_dipolar
end interface transmutation_visit_dipolar

contains

    pure subroutine visit_walls(overlap, delta_energies, visitable_walls, wall_pairs, &
        ij_components,  particles)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)

        integer :: i_partner

        do i_partner = size(particles), 1, -1
            call visitable_walls%visit(overlap, delta_energies(i_partner), &
                particles(i_partner)%position, wall_pairs(ij_components(i_partner))%potential)
            if (overlap) return
        end do
        delta_energies(1) = -delta_energies(1)
    end subroutine visit_walls

    subroutine visit_short(overlap, delta_energies, cells, ij_components, particles)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:, :)
        type(Cells_Wrapper), intent(in) :: cells
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)

        integer :: k_component, i_partner, i_exclude

        do i_partner = size(particles), 1, -1
            do k_component = 1, size(cells%visitable_cells, 1)
                i_exclude = i_exclude_particle(k_component, ij_components, particles)
                call cells%visitable_cells(k_component, ij_components(i_partner))%&
                    visit_energy(overlap, delta_energies(k_component, i_partner), &
                        particles(i_partner), visit_different, i_exclude)
                if (overlap) return
            end do
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
    end subroutine visit_short

    pure subroutine visit_field(delta_energies, external_field, particles)
        real(DP), intent(out) :: delta_energies(:)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Concrete_Particle), intent(in) :: particles(:)

        delta_energies = &
            [dipoles_field_visit_remove(external_field, particles(1)), &
            dipoles_field_visit_add(external_field, particles(2))]
    end subroutine visit_field

    pure subroutine visit_dipolar(delta_energies, delta_shared_energy, &
        dipolar_interactions_dynamic, ij_components, particles)
        real(DP), intent(out) :: delta_energies(:, :)
        real(DP), intent(out) :: delta_shared_energy
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)

        integer :: k_component, i_partner, i_exclude

        do i_partner = 1, size(particles)
            do k_component = 1, size(dipolar_interactions_dynamic%real_components, 1)
                i_exclude = i_exclude_particle(k_component, ij_components, particles)
                call dipolar_interactions_dynamic%&
                    real_components(k_component, ij_components(i_partner))%component%&
                        visit(delta_energies(k_component, i_partner), particles(i_partner), &
                            visit_different, i_exclude)
            end do
            delta_energies(ij_components(i_partner), i_partner) = &
                delta_energies(ij_components(i_partner), i_partner) - &
                dipolar_interactions_dynamic%self_components(ij_components(i_partner))%&
                component%meet(particles(i_partner)%dipole_moment)
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
        delta_shared_energy = &
            dipolar_interactions_dynamic%reci_visitor%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1)) + &
            dipolar_interactions_dynamic%surf_mixture%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1)%&
                    dipole_moment) - &
            dipolar_interactions_dynamic%dlc_visitor%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
    end subroutine visit_dipolar

    !> @note The i_actor <-> j_actor term is ignored.
    pure integer function i_exclude_particle(k_component, ij_components, particles) &
        result(i_exclude)
        integer, intent(in) :: k_component, ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)

        if (k_component == ij_components(1)) then
            i_exclude = particles(1)%i
        else if (k_component == ij_components(2)) then
            i_exclude = particles(2)%i
        else
            i_exclude = 0
        end if
    end function i_exclude_particle

end module procedures_swap_visitors
