module procedures_short_interactions_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_visit_condition, only: abstract_visit_condition, visit_lower, visit_all
use types_particle_wrapper, only: Concrete_Particle
use types_component_wrapper, only: Component_Wrapper
use classes_pair_potential, only: Pair_Potential_Wrapper, Pair_Potential_Line
use classes_walls_visitor, only: Abstract_Walls_Visitor
use classes_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor
use classes_visitable_cells, only: Abstract_Visitable_Cells
use types_real_wrapper, only: Real_Line

implicit none

private
public :: visit

interface visit
    module procedure :: visit_walls
    module procedure :: visit_short
    module procedure :: visit_cells_energies
    module procedure :: visit_cells_contacts
    module procedure :: visit_cells_min_distance
end interface visit

contains

    pure subroutine visit_walls(overlap, energies, components, walls_visitors, wall_pairs)
        logical, intent(out) :: overlap
        real(DP), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Walls_Visitor), intent(in) :: walls_visitors
        type(Pair_Potential_Wrapper), intent(in) :: wall_pairs(:)

        integer :: i_component

        overlap = .false.
        do i_component = 1, size(energies)
            call walls_visitors%visit(overlap, energies(i_component), &
                components(i_component)%positions, wall_pairs(i_component)%potential)
            if (overlap) return
        end do
    end subroutine visit_walls

    pure subroutine visit_short(overlap, energies, components, components_visitor, components_pairs)
        logical, intent(out) :: overlap
        type(Real_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Short_Pairs_Visitor), intent(in) :: components_visitor
        type(Pair_Potential_Line), intent(in) :: components_pairs(:)

        call visit_short_intra(overlap, energies, components, components_visitor, components_pairs)
        if (overlap) return
        call visit_short_inter(overlap, energies, components, components_visitor, components_pairs)
    end subroutine visit_short

    pure subroutine visit_short_intra(overlap, energies, components, components_visitor, &
        components_pairs)
        logical, intent(out) :: overlap
        type(Real_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Short_Pairs_Visitor), intent(in) :: components_visitor
        type(Pair_Potential_Line), intent(in) :: components_pairs(:)

        integer :: i_component

        overlap = .false.
        do i_component = 1, size(components)
            associate(energy_i => energies(i_component)%line(i_component), &
                positions_i => components(i_component)%positions, &
                potential_i => components_pairs(i_component)%line(i_component)%potential)
            call components_visitor%visit(overlap, energy_i, positions_i, potential_i)
            end associate
            if (overlap) return
        end do
    end subroutine visit_short_intra

    pure subroutine visit_short_inter(overlap, energies, components, components_visitor, &
        components_pairs)
        logical, intent(out) :: overlap
        type(Real_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Short_Pairs_Visitor), intent(in) :: components_visitor
        type(Pair_Potential_Line), intent(in) :: components_pairs(:)

        integer :: i_component, j_component

        overlap = .false.
        do j_component = 1, size(components)
            do i_component = 1, j_component - 1
                associate(energy_ij => energies(j_component)%line(i_component), &
                    positions_i => components(i_component)%positions, &
                    positions_j => components(j_component)%positions, &
                    potential_ij => components_pairs(j_component)%line(i_component)%potential)
                    call components_visitor%visit(overlap, energy_ij, positions_i, positions_j, &
                        potential_ij)
                end associate
                if (overlap) return
            end do
        end do
    end subroutine visit_short_inter

    subroutine visit_cells_energies(overlap, energies, components, visitable_cells)
        logical, intent(out) :: overlap
        type(Real_Line), intent(inout) :: energies(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells(:, :)

        real(DP) :: energy_ij, energy_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Particle) :: particle
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
                    call visitable_cells(i_component, j_component)%&
                        visit_energy(overlap, energy_j, particle, visit_condition, i_exclude)
                    if (overlap) return
                    energy_ij = energy_ij + energy_j
                end do
                energies(j_component)%line(i_component) = energy_ij
            end do
        end do
    end subroutine visit_cells_energies

    subroutine visit_cells_contacts(overlap, contacts, components, visitable_cells)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: contacts
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells(:, :)

        real(DP) :: conctacts_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Particle) :: particle
        procedure(abstract_visit_condition), pointer :: visit_condition => null()

        overlap = .false.
        contacts = 0._DP
        do j_component = 1, size(visitable_cells, 2)
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
                    call visitable_cells(i_component, j_component)%&
                        visit_contacts(overlap, conctacts_j, particle, visit_condition, i_exclude)
                    if (overlap) return
                    contacts = contacts + conctacts_j
                end do
            end do
        end do
    end subroutine visit_cells_contacts

    subroutine visit_cells_min_distance(overlap, min_distance_ratio, max_distance_ratio, &
        components, visitable_cells)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: min_distance_ratio
        real(DP), intent(in) :: max_distance_ratio
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Visitable_Cells), intent(in) :: visitable_cells(:, :)

        real(DP) :: min_distance_ratio_j
        integer :: i_component, j_component, i_particle, i_exclude
        logical :: same_component
        type(Concrete_Particle) :: particle
        procedure(abstract_visit_condition), pointer :: visit_condition => null()

        overlap = .false.
        min_distance_ratio = max_distance_ratio
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
                    call visitable_cells(i_component, j_component)%&
                        visit_min_distance(overlap, min_distance_ratio_j, particle, &
                            visit_condition, i_exclude)
                    if (overlap) return
                    if (min_distance_ratio_j < min_distance_ratio) &
                        min_distance_ratio = min_distance_ratio_j
                end do
            end do
        end do
    end subroutine visit_cells_min_distance

end module procedures_short_interactions_visitor
