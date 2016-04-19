module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use classes_number_to_string, only: Concrete_Number_to_String
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_external_field, only: Abstract_External_Field
use classes_walls_potential, only: Abstract_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use classes_pair_potential, only: Abstract_Pair_Potential
use classes_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: DES_Self_Component_Wrapper, &
    Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component
use types_line_observables, only: Concrete_Line_Observables
use types_observables_wrapper, only: Observables_Wrapper
use procedures_observables_factory, only: create_triangle_nodes, destroy_triangle_nodes
use procedures_triangle_observables, only: triangle_observables_init, &
    triangle_observables_add

implicit none

private
public :: plmc_visit

interface plmc_visit
    module procedure :: visit_all
    module procedure :: visit_field
    module procedure :: visit_walls
    module procedure :: visit_short
    module procedure :: visit_dipolar
end interface plmc_visit

contains

    subroutine visit_all(observables, environment, mixture, short_interactions, &
        dipolar_interactions)
        type(Observables_Wrapper), intent(inout) :: observables
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        call plmc_visit(observables%field_energies, environment%external_field, mixture%components)
        call plmc_visit(observables%walls_energies, mixture%components, short_interactions)
        call plmc_visit(observables%short_energies, mixture%components, short_interactions)
        call plmc_visit(observables%dipolar_energies, observables%dipolar_mixture_energy, mixture%&
            components, dipolar_interactions)
    end subroutine visit_all

    pure subroutine visit_field(field_energies, external_field, components)
        real(DP), intent(out) :: field_energies(:)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(components)
            field_energies(i_component) = dipoles_field_visit_component(external_field, &
                components(i_component)%positions, components(i_component)%dipolar_moments)
        end do
    end subroutine visit_field

    subroutine visit_walls(walls_energies, components, short_interactions)
        real(DP), intent(out) :: walls_energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        logical :: overlap
        integer :: i_component
        type(Concrete_Number_to_String) :: string

        do i_component = 1, size(components)
            associate(energy_i => walls_energies(i_component), &
                positions_i => components(i_component)%positions, &
                potential_i => short_interactions%wall_pairs(i_component)%potential)
                call short_interactions%walls_visitor%visit(overlap, energy_i, positions_i, &
                    potential_i)
            end associate
            if (overlap) then
                call error_exit("visit_walls: component "//string%get(i_component)//&
                    " overlaps.")
            end if
        end do
    end subroutine visit_walls

    subroutine visit_short(energies, components, short_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        call visit_short_intra(energies, components, short_interactions)
        call visit_short_inter(energies, components, short_interactions)
    end subroutine visit_short

    subroutine visit_short_intra(energies, components, short_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        logical :: overlap
        integer :: i_component
        type(Concrete_Number_to_String) :: string

        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%line(i_component), &
                positions_i => components(i_component)%positions, &
                potential_i => short_interactions%components_pairs(i_component)%line(i_component)%&
                    potential)
            call short_interactions%components_visitor%visit(overlap, energy_i, positions_i, &
                potential_i)
            end associate
            if (overlap) then
                call error_exit("visit_short_intra: component "//string%get(i_component)//&
                    " overlaps with itself.")
            end if
        end do
    end subroutine visit_short_intra

    subroutine visit_short_inter(energies, components, short_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        logical :: overlap
        integer :: i_component, j_component
        type(Concrete_Number_to_String) :: string

        do j_component = 1, size(components)
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%line(i_component), &
                    positions_i => components(i_component)%positions, &
                    positions_j => components(j_component)%positions, &
                    potential_ij => short_interactions%components_pairs(j_component)%&
                        line(i_component)%potential)
                    call short_interactions%components_visitor%visit(overlap, energy_ij, &
                        positions_i, positions_j, potential_ij)
                end associate
                if (overlap) then
                    call error_exit("visit_short_inter: components "//string%get(i_component)//&
                        " and "//string%get(j_component)//" overlap.")
                end if
            end do
        end do
    end subroutine visit_short_inter

    pure subroutine visit_dipolar(energies, dipolar_mixture_energy, components, &
        dipolar_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        real(DP), intent(out) :: dipolar_mixture_energy
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        type(Concrete_Line_Observables) :: real_energies(size(components))
        real(DP) :: self_energies(size(components))

        call triangle_observables_init(energies)

        call create_triangle_nodes(real_energies)
        call visit_des_real(real_energies, components, dipolar_interactions)
        call triangle_observables_add(energies, real_energies)
        call destroy_triangle_nodes(real_energies)

        dipolar_mixture_energy = dipolar_interactions%reci_visitor%visit() + dipolar_interactions%&
            surf_mixture%visit() - dipolar_interactions%dlc_visitor%visit()

        call visit_des_self(self_energies, dipolar_interactions%self_components)
        call triangle_observables_add(energies, -self_energies)
    end subroutine visit_dipolar

    pure subroutine visit_des_real(energies, components, dipolar_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        call visit_des_real_intra(energies, components, dipolar_interactions)
        call visit_des_real_inter(energies, components, dipolar_interactions)
    end subroutine visit_des_real

    pure subroutine visit_des_real_intra(energies, components, dipolar_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        integer :: i_component

        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%line(i_component), &
                    positions_i => components(i_component)%positions, &
                    dipolar_moments_i => components(i_component)%dipolar_moments, &
                    potential_i => dipolar_interactions%real_pairs(i_component)%line(i_component)%&
                        potential)
                call dipolar_interactions%real_visitor%visit(energy_i, positions_i, &
                    dipolar_moments_i, potential_i)
            end associate
        end do
    end subroutine visit_des_real_intra

    pure subroutine visit_des_real_inter(energies, components, dipolar_interactions)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        integer :: i_component, j_component

        do j_component = 1, size(components)
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%line(i_component), &
                    positions_i => components(i_component)%positions, &
                    dipolar_moments_i => components(i_component)%dipolar_moments, &
                    positions_j => components(j_component)%positions, &
                    dipolar_moments_j => components(j_component)%dipolar_moments, &
                    potential_ij => dipolar_interactions%real_pairs(j_component)%line(i_component)%&
                        potential)
                    call dipolar_interactions%real_visitor%visit(energy_ij, positions_i, &
                        dipolar_moments_i, positions_j, dipolar_moments_j, potential_ij)
                end associate
            end do
        end do
    end subroutine visit_des_real_inter

    pure subroutine visit_des_self(energies, self_components)
        real(DP), intent(out) :: energies(:)
        type(DES_Self_Component_Wrapper), intent(in) :: self_components(:)

        integer :: i_component

        do i_component = 1, size(energies)
            energies(i_component) = self_components(i_component)%component%visit()
        end do
    end subroutine visit_des_self

end module procedures_plmc_visit
