module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_external_field, only: Abstract_External_Field
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_condition_lower => lower
use types_des_self_component_wrapper, only: DES_Self_Component_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_component => visit_component
use types_reals_line, only: Reals_Line
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use procedures_observables_factory_micro, only: create_triangle_nodes, destroy_triangle_nodes
use procedures_triangle_observables, only: triangle_observables_init, &
    triangle_observables_add

implicit none

private
public :: plmc_visit_set, plmc_visit, visit_short_full, visit_short_cells

interface plmc_visit_set
    module procedure :: set_visit
end interface plmc_visit_set

interface plmc_visit
    module procedure :: visit_generating, visit_exploring
    module procedure :: visit_field
    module procedure :: visit_walls
    module procedure :: visit_dipolar
end interface plmc_visit

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

    subroutine visit_generating(observables, physical_model)
        type(Generating_Observables_Wrapper), intent(inout) :: observables
        type(Physical_Model_Wrapper), intent(in) :: physical_model

        call plmc_visit(observables%field_energies, physical_model%environment%external_field, &
            physical_model%mixture%components)
        call plmc_visit(observables%walls_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        call visit_short_full(observables%short_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        call plmc_visit(observables%dipolar_energies, observables%dipolar_mixture_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions)
    end subroutine visit_generating

    subroutine visit_exploring(observables, physical_model, visit)
        type(Exploring_Observables_Wrapper), intent(inout) :: observables
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        logical, intent(in) :: visit

        if (.not.visit) return
        call plmc_visit(observables%field_energies, physical_model%environment%external_field, &
            physical_model%mixture%components)
        call plmc_visit(observables%walls_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        call visit_short_cells(observables%short_energies, physical_model%mixture%components, &
            physical_model%short_interactions)
        call plmc_visit(observables%dipolar_energies, observables%dipolar_mixture_energy, &
            physical_model%mixture%components, physical_model%dipolar_interactions)
    end subroutine visit_exploring

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
                call error_exit("procedures_plmc_visit: visit_walls: component "//string%&
                    get(i_component)//" overlaps.")
            end if
        end do
    end subroutine visit_walls

    subroutine visit_short_full(energies, components, short_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        call visit_short_intra(energies, components, short_interactions)
        call visit_short_inter(energies, components, short_interactions)
    end subroutine visit_short_full

    subroutine visit_short_intra(energies, components, short_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        integer :: i_component
        logical :: overlap
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
                call error_exit("procedures_plmc_visit: visit_short_intra: component "//string%&
                    get(i_component)//" overlaps with itself.")
            end if
        end do
    end subroutine visit_short_intra

    subroutine visit_short_inter(energies, components, short_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        integer :: i_component, j_component
        logical :: overlap
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
                    call error_exit("procedures_plmc_visit: visit_short_inter: components "//&
                        string%get(i_component)//" and "//string%get(j_component)//" overlap.")
                end if
            end do
        end do
    end subroutine visit_short_inter

    subroutine visit_short_cells(energies, components, short_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions

        real(DP) :: energy_ij, energy_j
        integer :: j_component, i_component, num_particles_max, i_particle, i_exclude
        logical :: overlap
        type(Concrete_Temporary_Particle) :: particle
        type(Concrete_Number_to_String) :: string

        do j_component = 1, size(short_interactions%visitable_cells, 2)
            do i_component = 1, j_component
                num_particles_max = max(components(j_component)%positions%get_num(), &
                    components(i_component)%positions%get_num())
                energy_ij = 0._DP
                do i_particle = 1, components(j_component)%positions%get_num()
                    particle%i = i_particle
                    particle%position = components(j_component)%positions%get(particle%i)
                    i_exclude = merge(particle%i, num_particles_max + 1, i_component == j_component)
                    call short_interactions%visitable_cells(i_component, j_component)%&
                        visit(overlap, energy_j, particle, visit_condition_lower, i_exclude)
                    if (overlap) then
                        call error_exit("procedures_plmc_visit: visit_short_cells: components "//&
                            string%get(i_component)//" and "//string%get(j_component)//" overlap.")
                    end if
                    energy_ij = energy_ij + energy_j
                end do
                energies(j_component)%line(i_component) = energy_ij
            end do
        end do
    end subroutine visit_short_cells

    pure subroutine visit_dipolar(energies, mixture_energy, components, &
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

    pure subroutine visit_des_real(energies, components, dipolar_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        call visit_des_real_intra(energies, components, dipolar_interactions)
        call visit_des_real_inter(energies, components, dipolar_interactions)
    end subroutine visit_des_real

    pure subroutine visit_des_real_intra(energies, components, dipolar_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        integer :: i_component

        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%line(i_component), &
                    positions_i => components(i_component)%positions, &
                    dipolar_moments_i => components(i_component)%dipolar_moments, &
                    pair => dipolar_interactions%real_pair)
                call dipolar_interactions%real_visitor%visit(energy_i, positions_i, &
                    dipolar_moments_i, pair)
            end associate
        end do
    end subroutine visit_des_real_intra

    pure subroutine visit_des_real_inter(energies, components, dipolar_interactions)
        type(Reals_Line), intent(inout) :: energies(:)
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
                    pair => dipolar_interactions%real_pair)
                    call dipolar_interactions%real_visitor%visit(energy_ij, positions_i, &
                        dipolar_moments_i, positions_j, dipolar_moments_j, pair)
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
