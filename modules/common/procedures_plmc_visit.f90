module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use class_number_to_string, only: Concrete_Number_to_String
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_walls_potential, only: Abstract_Walls_Potential
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use class_pair_potential, only: Abstract_Pair_Potential
use class_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use procedures_ewald_reci_visit, only: ewald_reci_visit
use procedures_ewald_self_visit, only: ewald_self_visit
use types_long_interactions_wrapper, only: Ewald_Self_Wrapper, Long_Interactions_Wrapper
use types_observables_wrapper, only: Concrete_Components_Energies, Observables_Wrapper
use procedures_observables_factory, only: create_components_energies_nodes, &
    destroy_components_energies_nodes
use procedures_components_energies, only: Concrete_Components_Energies_init, &
    Concrete_Components_Energies_add

implicit none

private
public :: plmc_visit

interface plmc_visit
    module procedure :: visit_all
    module procedure :: visit_walls
    module procedure :: visit_short
    module procedure :: visit_long
end interface plmc_visit

contains

    subroutine visit_all(observables, mixture, short_interactions, long_interactions)
        type(Observables_Wrapper), intent(inout) :: observables
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        call plmc_visit(observables%walls_energies, mixture%components, short_interactions)
        call plmc_visit(observables%short_energies, mixture%components, short_interactions)
        call plmc_visit(observables%long_energies, mixture%components, long_interactions)
    end subroutine visit_all

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
                pair_i => short_interactions%wall_pairs(i_component)%pair_potential)
                call short_interactions%walls_visitor%visit(overlap, energy_i, positions_i, pair_i)
            end associate
            if (overlap) then
                call error_exit("visit_walls: component "//string%get(i_component)//&
                    " overlaps.")
            end if
        end do
    end subroutine visit_walls

    subroutine visit_short(energies, components, short_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        call visit_short_intra(energies, components, short_interactions)
        call visit_short_inter(energies, components, short_interactions)
    end subroutine visit_short

    subroutine visit_short_intra(energies, components, short_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        logical :: overlap
        integer :: i_component
        type(Concrete_Number_to_String) :: string

        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%with_components(i_component), &
                positions_i => components(i_component)%positions, &
                pair_i => short_interactions%components_pairs(i_component)%&
                    with_components(i_component)%pair_potential)
            call short_interactions%components_visitor%visit(overlap, energy_i, positions_i, pair_i)
            end associate
            if (overlap) then
                call error_exit("visit_short_intra: component "//string%get(i_component)//&
                    " overlaps with itself.")
            end if
        end do
    end subroutine visit_short_intra

    subroutine visit_short_inter(energies, components, short_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        logical :: overlap
        integer :: j_component, i_component
        type(Concrete_Number_to_String) :: string

        do j_component = 1, size(components)
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%with_components(i_component), &
                    positions_i => components(i_component)%positions, &
                    positions_j => components(j_component)%positions, &
                    pair_ij => short_interactions%components_pairs(j_component)%&
                        with_components(i_component)%pair_potential)
                    call short_interactions%components_visitor%visit(overlap, energy_ij, &
                        positions_i, positions_j, pair_ij)
                end associate
                if (overlap) then
                    call error_exit("visit_short_inter: components "//string%get(i_component)//&
                        " and "//string%get(j_component)//" overlap.")
                end if
            end do
        end do
    end subroutine visit_short_inter

    pure subroutine visit_long(energies, components, long_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        type(Concrete_Components_Energies) :: real_energies(size(components))
        type(Concrete_Components_Energies) :: reci_energies(size(components))
        real(DP) :: self_energies(size(components))

        call Concrete_Components_Energies_init(energies)

        call create_components_energies_nodes(real_energies)
        call visit_long_real(real_energies, components, long_interactions)
        call Concrete_Components_Energies_add(energies, real_energies)
        call destroy_components_energies_nodes(real_energies)

        call create_components_energies_nodes(reci_energies)
        call visit_long_reci(reci_energies, long_interactions)
        call Concrete_Components_Energies_add(energies, reci_energies)
        call destroy_components_energies_nodes(reci_energies)

        call visit_long_self(self_energies, components, long_interactions%selves)
        call Concrete_Components_Energies_add(energies, -self_energies)
    end subroutine visit_long

    pure subroutine visit_long_real(energies, components, long_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        call visit_long_real_intra(energies, components, long_interactions)
        call visit_long_real_inter(energies, components, long_interactions)
    end subroutine visit_long_real

    pure subroutine visit_long_real_intra(energies, components, long_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        integer :: i_component

        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%with_components(i_component), &
                    positions_i => components(i_component)%positions, &
                    dipolar_moment_i => components(i_component)%dipolar_moments, &
                    pair_i => long_interactions%real_pairs(i_component)%&
                        with_components(i_component)%real_pair)
                call long_interactions%real_visitor%visit(energy_i, positions_i, dipolar_moment_i, &
                    pair_i)
            end associate
        end do
    end subroutine visit_long_real_intra

    pure subroutine visit_long_real_inter(energies, components, long_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        integer :: j_component, i_component

        do j_component = 1, size(components)
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%with_components(i_component), &
                    positions_i => components(i_component)%positions, &
                    dipolar_moments_i => components(i_component)%dipolar_moments, &
                    positions_j => components(j_component)%positions, &
                    dipolar_moments_j => components(j_component)%dipolar_moments, &
                    pair_ij => long_interactions%real_pairs(j_component)%&
                        with_components(i_component)%real_pair)
                    call long_interactions%real_visitor%visit(energy_ij, positions_i, &
                        dipolar_moments_i, positions_j, dipolar_moments_j, pair_ij)
                end associate
            end do
        end do
    end subroutine visit_long_real_inter

    pure subroutine visit_long_reci(energies, long_interactions)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        integer :: j_component, i_component

        do j_component = 1, size(energies)
            associate(energy_j => energies(j_component)%with_components(j_component), &
                structure_j => long_interactions%reci_components(j_component)%reci_component)
                energy_j = ewald_reci_visit(long_interactions%reci_weight, structure_j)
            end associate
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%with_components(i_component), &
                    structure_i => long_interactions%reci_components(i_component)%reci_component, &
                    structure_j => long_interactions%reci_components(j_component)%reci_component)
                    energy_ij = ewald_reci_visit(long_interactions%reci_weight, structure_i, &
                        structure_j)
                end associate
            end do
        end do
    end subroutine visit_long_reci

    pure subroutine visit_long_self(energies, components, ewald_selves)
        real(DP), intent(out) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Ewald_Self_Wrapper), intent(in) :: ewald_selves(:)

        integer :: i_component

        do i_component = 1, size(energies)
            energies(i_component) = ewald_self_visit(components(i_component)%dipolar_moments, &
                ewald_selves(i_component)%self)
        end do
    end subroutine visit_long_self

end module procedures_plmc_visit
