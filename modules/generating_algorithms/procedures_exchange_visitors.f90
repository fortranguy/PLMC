module procedures_exchange_visitors

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_visit_condition, only: visit_different
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_cells_wrapper, only: Cells_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper

implicit none

private
public :: visit_add, visit_remove

interface visit_add
    module procedure :: visit_add_short
    module procedure :: visit_add_dipolar
end interface visit_add

interface visit_remove
    module procedure :: visit_remove_dipolar
    module procedure :: visit_remove_short
end interface visit_remove

contains

    subroutine visit_add_short(overlap, delta_energies, i_component, particle, cells)
        logical, intent(out) :: overlap
        real(DP), intent(inout) :: delta_energies(:)
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        type(Cells_Wrapper), intent(in) :: cells

        integer :: j_component, i_exclude

        do j_component = 1, size(cells%visitable_cells, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call cells%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, delta_energies(j_component), particle, visit_different, &
                    i_exclude)
            if (overlap) return
        end do
    end subroutine visit_add_short

    subroutine visit_remove_short(overlap, delta_energies, i_component, particle, cells)
        logical, intent(out) :: overlap
        real(DP), intent(inout) :: delta_energies(:)
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        type(Cells_Wrapper), intent(in) :: cells

        integer :: j_component, i_exclude

        do j_component = 1, size(cells%visitable_cells, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call cells%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, delta_energies(j_component), particle, visit_different, &
                    i_exclude)
        end do
        delta_energies = -delta_energies
    end subroutine visit_remove_short

    subroutine visit_add_dipolar(delta_energies, delta_shared_energy, i_component, particle, &
        dipolar_interactions_dynamic)
        real(DP), intent(inout) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic

        integer :: j_component, i_exclude

        do j_component = 1, size(dipolar_interactions_dynamic%real_components, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call dipolar_interactions_dynamic%real_components(j_component, i_component)%component%&
                visit(delta_energies(j_component), particle, visit_different, i_exclude)
        end do
        delta_energies(i_component) = delta_energies(i_component) - dipolar_interactions_dynamic%&
            self_components(i_component)%component%meet(particle%dipole_moment)

        delta_shared_energy = &
            dipolar_interactions_dynamic%reci_visitor%visit_add(i_component, particle) + &
            dipolar_interactions_dynamic%surf_mixture%&
                visit_add(i_component, particle%dipole_moment) - &
            dipolar_interactions_dynamic%dlc_visitor%visit_add(i_component, particle)
    end subroutine visit_add_dipolar

    subroutine visit_remove_dipolar(delta_energies, delta_shared_energy, i_component, particle, &
        dipolar_interactions_dynamic)
        real(DP), intent(inout) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic

        integer :: j_component, i_exclude

        do j_component = 1, size(dipolar_interactions_dynamic%real_components, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call dipolar_interactions_dynamic%real_components(j_component, i_component)%component%&
                visit(delta_energies(j_component), particle, visit_different, i_exclude)
        end do
        delta_energies(i_component) = delta_energies(i_component) - dipolar_interactions_dynamic%&
            self_components(i_component)%component%meet(particle%dipole_moment)
        delta_energies = -delta_energies

        delta_shared_energy = &
            dipolar_interactions_dynamic%reci_visitor%visit_remove(i_component, particle) + &
            dipolar_interactions_dynamic%surf_mixture%&
                visit_remove(i_component, particle%dipole_moment) - &
            dipolar_interactions_dynamic%dlc_visitor%visit_remove(i_component, particle)
    end subroutine visit_remove_dipolar

end module procedures_exchange_visitors
