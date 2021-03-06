module procedures_dipolar_interactions_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_visit_condition, only: abstract_visit_condition, visit_lower, visit_all
use classes_external_field, only: Abstract_External_Field
use types_particle_wrapper, only: Concrete_Particle
use types_component_wrapper, only: Component_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component => visit_component
use classes_des_real_component, only: DES_Real_Component_Wrapper
use classes_des_self_component, only: DES_Self_Component_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_real_wrapper, only: Real_Line
use procedures_reals_factory, only: reals_create => create, reals_destroy => destroy
use procedures_triangle_observables, only: triangle_observables_init, triangle_observables_add

implicit none

private
public :: visit

interface visit
    module procedure :: visit_field
    module procedure :: visit_dipolar
    module procedure :: visit_des_real
end interface visit

contains

    pure subroutine visit_field(energies, external_field, components)
        real(DP), intent(inout) :: energies(:)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(energies)
            energies(i_component) = dipoles_field_visit_component(external_field, &
                components(i_component)%positions, components(i_component)%dipole_moments)
        end do
    end subroutine visit_field

    subroutine visit_dipolar(energies, shared_energy, components, dipolar_interactions_dynamic)
        type(Real_Line), intent(inout) :: energies(:)
        real(DP), intent(out) :: shared_energy
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic

        type(Real_Line), allocatable :: real_energies(:)
        real(DP) :: self_energies(size(components))

        call triangle_observables_init(energies)

        call reals_create(real_energies, size(components))
        call visit_des_real(real_energies, components, dipolar_interactions_dynamic%real_components)
        call triangle_observables_add(energies, real_energies)
        call reals_destroy(real_energies)

        call visit_des_self(self_energies, dipolar_interactions_dynamic%self_components)
        call triangle_observables_add(energies, -self_energies)

        shared_energy = dipolar_interactions_dynamic%reci_visitor%visit() + &
            dipolar_interactions_dynamic%surf_mixture%visit() - &
            dipolar_interactions_dynamic%dlc_visitor%visit()
    end subroutine visit_dipolar

    subroutine visit_des_real(energies, components, real_components)
        type(Real_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(DES_Real_Component_Wrapper), intent(in) :: real_components(:, :)

        real(DP) :: energy_ij, energy_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Particle) :: particle
        procedure(abstract_visit_condition), pointer :: visit_condition => null()

        do j_component = 1, size(energies)
            do i_component = 1, j_component
                same_component = i_component == j_component
                if (same_component) then
                    visit_condition => visit_lower
                else
                    visit_condition => visit_all
                end if
                energy_ij = 0._DP
                do i_particle = 1, components(j_component)%dipole_moments%get_num()
                    particle%i = i_particle
                    particle%position = components(j_component)%positions%get(particle%i)
                    particle%dipole_moment = components(j_component)%dipole_moments%get(particle%i)
                    i_exclude = merge(particle%i, 0, same_component)
                    call real_components(i_component, j_component)%component%&
                        visit(energy_j, particle, visit_condition, i_exclude)
                    energy_ij = energy_ij + energy_j
                end do
                energies(j_component)%line(i_component) = energy_ij
            end do
        end do
    end subroutine visit_des_real

    pure subroutine visit_des_self(energies, self_components)
        real(DP), intent(inout) :: energies(:)
        type(DES_Self_Component_Wrapper), intent(in) :: self_components(:)

        integer :: i_component

        do i_component = 1, size(energies)
            energies(i_component) = self_components(i_component)%component%visit()
        end do
    end subroutine visit_des_self

end module procedures_dipolar_interactions_visitor
