module procedures_plmc_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use class_walls_potential, only: Abstract_Walls_Potential
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Mixture_Wrapper_Old
use class_pair_potential, only: Abstract_Pair_Potential
use class_component_potential, only: Abstract_Component_Potential
use types_short_potentials_wrapper, only: Mixture_Short_Potentials_Wrapper
use class_ewald_real_component, only: Abstract_Ewald_Real_Component
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use types_observables_wrapper, only: Observables_Wrapper

implicit none

private
public :: plmc_visit

interface plmc_visit
    module procedure :: visit_mixture
    module procedure :: visit_mixture_short
    module procedure :: visit_component_walls
    module procedure :: visit_component_short
    module procedure :: visit_mixture_ewald
    module procedure :: visit_component_ewald_real
end interface plmc_visit

contains

    subroutine visit_mixture(observables, walls_potential, short_potentials, ewalds, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        type(Mixture_Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Mixture_Wrapper_Old), intent(in) :: mixture

        call plmc_visit(observables, walls_potential, short_potentials, mixture)
        call plmc_visit(observables, ewalds, mixture)
    end subroutine visit_mixture

    subroutine visit_mixture_short(observables, walls_potential, short_potentials, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        type(Mixture_Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Wrapper_Old), intent(in) :: mixture

        logical :: overlap
        real(DP) :: inter_energy_1, inter_energy_2

        call plmc_visit(overlap, observables%intras(1)%component_energy%walls, walls_potential, &
            mixture%components(1)%positions, short_potentials%intras(1)%wall_pair)
        if (overlap) call error_exit("walls - components(1) overlap")
        call plmc_visit(overlap, observables%intras(2)%component_energy%walls, walls_potential, &
            mixture%components(2)%positions, short_potentials%intras(2)%wall_pair)
        if (overlap) call error_exit("walls - components(2) overlap")

        call plmc_visit(overlap, observables%intras(1)%component_energy%short, &
            short_potentials%intras(1)%component, mixture%components(1)%positions, &
            short_potentials%intras(1)%pair, same_type=.true.)
        if (overlap) call error_exit("intra 1 overlap")
        call plmc_visit(overlap, observables%intras(2)%component_energy%short, &
            short_potentials%intras(2)%component, mixture%components(2)%positions, &
            short_potentials%intras(2)%pair, same_type=.true.)
        if (overlap) call error_exit("intra 2 overlap")

        call plmc_visit(overlap, inter_energy_1, short_potentials%intras(1)%component, &
            mixture%components(2)%positions, short_potentials%inter_pair, same_type=.false.)
        if (overlap) call error_exit("inter intras(1) overlap")
        call plmc_visit(overlap, inter_energy_2, short_potentials%intras(2)%component, &
            mixture%components(1)%positions, short_potentials%inter_pair, same_type=.false.)
        if (overlap) call error_exit("inter intras(2) overlap")
        observables%inter_energy%short = inter_energy_1 + inter_energy_2
    end subroutine visit_mixture_short

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

    pure subroutine visit_component_short(overlap, energy, potential, positions, pair, same_type)
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Potential), intent(in) :: potential
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair
        logical, intent(in) :: same_type

        type(Concrete_Temporary_Particle) :: particle
        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, positions%get_num()
            particle%i = i_particle
            particle%position = positions%get(particle%i)
            call potential%visit(overlap, energy_i, particle, pair, same_type)
            if (overlap) return
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit_component_short

    pure subroutine visit_mixture_ewald(observables, ewalds, mixture)
        type(Observables_Wrapper), intent(inout) :: observables
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Mixture_Wrapper_Old), intent(in) :: mixture

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
            particle%i = i_particle
            particle%position = positions%get(particle%i)
            particle%dipolar_moment = dipolar_moments%get(particle%i)
            call potential%visit(energy_i, particle, same_type)
            energy = energy + energy_i
        end do
        energy = energy / 2._DP
    end subroutine visit_component_ewald_real

end module procedures_plmc_visit
