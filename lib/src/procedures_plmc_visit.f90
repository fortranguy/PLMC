module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use class_number_to_string, only: Concrete_Number_to_String
use class_walls_potential, only: Abstract_Walls_Potential
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_mixture_wrapper, only: Mixture_Wrapper
use class_pair_potential, only: Abstract_Pair_Potential
use class_short_potential_visitor, only: Abstract_Short_Potential_Visitor
use types_short_potentials_wrapper, only: Short_Potentials_Wrapper
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use types_ewalds_wrapper, only: Mixture_Ewald_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: plmc_visit

interface plmc_visit
    module procedure :: visit_mixture
    module procedure :: visit_mixture_short
    module procedure :: visit_component_walls
    module procedure :: visit_mixture_ewald
    module procedure :: visit_component_ewald_real
end interface plmc_visit

contains

    subroutine visit_mixture(observables, walls_potential, short_potentials, ewalds, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        type(Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Mixture_Wrapper), intent(in) :: mixture

        call plmc_visit(observables, short_potentials, mixture)
        call plmc_visit(observables, ewalds, mixture)
    end subroutine visit_mixture

    subroutine visit_mixture_short(observables, short_potentials, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        type(Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Wrapper), intent(in) :: mixture

        logical :: overlap
        integer :: j_component, i_component
        type(Concrete_Number_to_String) :: string

        do i_component = 1, size(mixture%components)
            associate(energy_i => observables%short_energies(i_component)%&
                with_components(i_component), &
                positions_i => mixture%components(i_component)%positions, &
                pair_i => short_potentials%inter_pairs(i_component)%with_components(i_component)%&
                    pair_potential)
            call short_potentials%visitor%visit(overlap, energy_i, positions_i, pair_i)
            end associate
            if (overlap) then
                call error_exit("visit_mixture_short: component "//string%get(i_component)//&
                    " overlaps with itself.")
            end if
        end do
        do j_component = 1, size(mixture%components)
            do i_component = 1, j_component - 1
                associate(energy_ij => observables%short_energies(j_component)%&
                    with_components(i_component), &
                    positions_i => mixture%components(i_component)%positions, &
                    positions_j => mixture%components(j_component)%positions, &
                    pair_ij => short_potentials%inter_pairs(j_component)%&
                        with_components(i_component)%pair_potential)
                    call short_potentials%visitor%visit(overlap, energy_ij, positions_i, &
                        positions_j, pair_ij)
                end associate
                if (overlap) then
                    call error_exit("visit_mixture_short: components "//string%get(i_component)//&
                        " and "//string%get(j_component)//" overlap.")
                end if
            end do
        end do
    end subroutine visit_mixture_short

        !call plmc_visit(overlap, observables%intras(1)%component_energy%walls, walls_potential, &
        !    mixture%components(1)%positions, short_potentials%intras(1)%wall_pair)
        !if (overlap) call error_exit("walls - components(1) overlap")

    pure subroutine visit_component_walls(overlap, energy, potential, positions, pair)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Walls_Potential), intent(in) :: potential
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair

        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, positions%get_num()
            call potential%visit(overlap, energy_i, positions%get(i_particle), pair)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine visit_component_walls

    pure subroutine visit_mixture_ewald(observables, ewalds, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Mixture_Wrapper), intent(in) :: mixture

        real(DP) :: real_energy_1, real_energy_2
        real(DP) :: inter_energy_1, inter_energy_2

        call plmc_visit(real_energy_1, ewalds%intras(1)%real_component, &
            mixture%components(1)%positions, mixture%components(1)%dipolar_moments, &
            same_type=.true.)
        call plmc_visit(real_energy_2, ewalds%intras(2)%real_component, &
            mixture%components(2)%positions, mixture%components(2)%dipolar_moments, &
            same_type=.true.)

        call plmc_visit(inter_energy_1, ewalds%inters(1)%real_component, &
            mixture%components(2)%positions, mixture%components(2)%dipolar_moments, &
            same_type=.false.)
        call plmc_visit(inter_energy_2, ewalds%inters(2)%real_component, &
            mixture%components(1)%positions, mixture%components(1)%dipolar_moments, &
            same_type=.false.)

        observables%intras(1)%component_energy%long = real_energy_1
        observables%intras(2)%component_energy%long = real_energy_2
        observables%inter_energy%long = inter_energy_1 + inter_energy_2
    end subroutine visit_mixture_ewald

    pure subroutine visit_component_ewald_real(energy, potential, positions, dipolar_moments, &
        same_type)
        real(DP), intent(out) :: energy
        class(Abstract_Ewald_Real_Component), intent(in) :: potential
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        logical, intent(in) :: same_type

        type(Concrete_Temporary_Particle) :: particle
        real(DP) :: energy_i
        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, positions%get_num()
            particle%position = positions%get(particle%i)
            particle%dipolar_moment = dipolar_moments%get(particle%i)
            !call potential%visit(energy_i, particle, i_particle) pair needed
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit_component_ewald_real

end module procedures_plmc_visit
