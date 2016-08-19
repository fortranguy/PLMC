module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_external_field, only: Abstract_External_Field
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: abstract_visit_condition, visit_lower, visit_all
use types_des_self_component_wrapper, only: DES_Self_Component_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component => visit_component
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_reals_line, only: Reals_Line
use types_observables_energies, only: Concrete_Energies
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use procedures_observables_factory_micro, only: create_triangle_nodes, destroy_triangle_nodes
use procedures_triangle_observables, only: triangle_observables_init, &
    triangle_observables_add

implicit none

private
public :: plmc_visit_set, plmc_visit, visit_walls, visit_short, visit_field, visit_dipolar

interface plmc_visit_set
    module procedure :: set_visit
end interface plmc_visit_set

interface plmc_visit
    module procedure :: visit_generating, visit_exploring
end interface plmc_visit

interface visit_short
    module procedure :: visit_short_energies
    module procedure :: visit_short_contacts
end interface visit_short

contains

    subroutine set_visit(visit, exploring_data, prefix)
        logical, intent(out) :: visit
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"visit"
        call exploring_data%get(data_field, visit, data_found)
        call check_data_found(data_field, data_found)
    end subroutine set_visit

    subroutine visit_generating(energies, physical_model)
        type(Concrete_Energies), intent(inout) :: energies
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        logical :: overlap

        call visit_walls(overlap, energies%walls_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visit: visit_generating: visit_walls: "&
            //"overlap.")
        call visit_short_full(overlap, energies%short_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visit: visit_generating: visit_short_full: "&
            //"overlap.")
        call visit_field(energies%field_energies, physical_model%environment%external_field, &
            physical_model%mixture%components)
        call visit_dipolar(energies%dipolar_energies, energies%dipolar_mixture_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions)
    end subroutine visit_generating

    subroutine visit_exploring(energies, physical_model, visit)
        type(Concrete_Energies), intent(inout) :: energies
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, intent(in) :: visit

        logical :: overlap

        if (.not.visit) return
        call visit_walls(overlap, energies%walls_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visit: visit_exploring: visit_walls: "&
            //"overlap.")
        call visit_short(overlap, energies%short_energies, physical_model%mixture%&
            components,physical_model%short_interactions)
        if (overlap) call error_exit("procedures_plmc_visit: visit_exploring: "//&
            "visit_short: overlap.")
        call visit_field(energies%field_energies, physical_model%environment%external_field, &
            physical_model%mixture%components)
        call visit_dipolar(energies%dipolar_energies, energies%dipolar_mixture_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions)
    end subroutine visit_exploring

    pure subroutine visit_walls(overlap, walls_energies, components, short_interactions)
        logical, intent(out) :: overlap
        real(DP), intent(inout) :: walls_energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        integer :: i_component

        overlap = .false.
        do i_component = 1, size(walls_energies)
            call short_interactions%walls_visitor%visit(overlap, walls_energies(i_component), &
                components(i_component)%positions, &
                short_interactions%wall_pairs(i_component)%potential)
            if (overlap) return
        end do
    end subroutine visit_walls

    pure subroutine visit_short_full(overlap, energies, components, short_interactions)
        logical, intent(out) :: overlap
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        call visit_short_intra(overlap, energies, components, short_interactions)
        if (overlap) return
        call visit_short_inter(overlap, energies, components, short_interactions)
    end subroutine visit_short_full

    pure subroutine visit_short_intra(overlap, energies, components, short_interactions)
        logical, intent(out) :: overlap
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        integer :: i_component

        overlap = .false.
        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%line(i_component), &
                positions_i => components(i_component)%positions, &
                potential_i => short_interactions%components_pairs(i_component)%line(i_component)%&
                    potential)
            call short_interactions%components_visitor%visit(overlap, energy_i, positions_i, &
                potential_i)
            end associate
            if (overlap) return
        end do
    end subroutine visit_short_intra

    pure subroutine visit_short_inter(overlap, energies, components, short_interactions)
        logical, intent(out) :: overlap
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        integer :: i_component, j_component

        overlap = .false.
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
                if (overlap) return
            end do
        end do
    end subroutine visit_short_inter

    subroutine visit_short_energies(overlap, energies, components, short_interactions)
        logical, intent(out) :: overlap
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        real(DP) :: energy_ij, energy_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Temporary_Particle) :: particle
        procedure(abstract_visit_condition), pointer :: visit_condition => null()

        overlap = .false.
        do j_component = 1, size(energies)
            do i_component = 1, j_component
                same_component = i_component == j_component
                if (same_component) then
                    visit_condition => visit_lower
                else
                    visit_condition => visit_all
                end if
                energy_ij = 0._DP
                do i_particle = 1, components(j_component)%positions%get_num()
                    particle%i = i_particle
                    particle%position = components(j_component)%positions%get(particle%i)
                    i_exclude = merge(particle%i, 0, same_component)
                    call short_interactions%visitable_cells(i_component, j_component)%&
                        visit_energy(overlap, energy_j, particle, visit_condition, i_exclude)
                    if (overlap) return
                    energy_ij = energy_ij + energy_j
                end do
                energies(j_component)%line(i_component) = energy_ij
            end do
        end do
    end subroutine visit_short_energies

    subroutine visit_short_contacts(overlap, contacts, components, short_interactions)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        real(DP) :: conctacts_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Temporary_Particle) :: particle
        procedure(abstract_visit_condition), pointer :: visit_condition => null()

        overlap = .false.
        contacts = 0._DP
        do j_component = 1, size(components)
            do i_component = 1, j_component
                same_component = i_component == j_component
                if (same_component) then
                    visit_condition => visit_lower
                else
                    visit_condition => visit_all
                end if
                do i_particle = 1, components(j_component)%positions%get_num()
                    particle%i = i_particle
                    particle%position = components(j_component)%positions%get(particle%i)
                    i_exclude = merge(particle%i, 0, same_component)
                    call short_interactions%visitable_cells(i_component, j_component)%&
                        visit_contacts(overlap, conctacts_j, particle, visit_condition, i_exclude)
                    if (overlap) return
                    contacts = contacts + conctacts_j
                end do
            end do
        end do
    end subroutine visit_short_contacts

    pure subroutine visit_field(field_energies, external_field, components)
        real(DP), intent(inout) :: field_energies(:)
        class(Abstract_External_Field), intent(in) :: external_field
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(field_energies)
            field_energies(i_component) = dipoles_field_visit_component(external_field, &
                components(i_component)%positions, components(i_component)%dipole_moments)
        end do
    end subroutine visit_field

    subroutine visit_dipolar(energies, mixture_energy, components, &
        dipolar_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        real(DP), intent(out) :: mixture_energy
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        type(Reals_Line) :: real_energies(size(components))
        real(DP) :: self_energies(size(components))

        call triangle_observables_init(energies)

        call create_triangle_nodes(real_energies)
        call visit_des_real(real_energies, components, dipolar_interactions)
        call triangle_observables_add(energies, real_energies)
        call destroy_triangle_nodes(real_energies)

        mixture_energy = dipolar_interactions%reci_visitor%visit() + dipolar_interactions%&
            surf_mixture%visit() - dipolar_interactions%dlc_visitor%visit()

        call visit_des_self(self_energies, dipolar_interactions%self_components)
        call triangle_observables_add(energies, -self_energies)
    end subroutine visit_dipolar

    subroutine visit_des_real(energies, components, dipolar_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        real(DP) :: energy_ij, energy_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Temporary_Particle) :: particle
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
                    call dipolar_interactions%real_components(i_component, j_component)%component%&
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

end module procedures_plmc_visit
