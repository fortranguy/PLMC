module procedures_observables_energies_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line
use types_observables_energies, only: Concrete_Single_Energies, Concrete_Double_Energies, &
    Concrete_Energies
use procedures_observables_factory_micro, only: create_reals, create_triangle_reals, destroy_reals,&
    destroy_triangle_reals

implicit none

private
public :: create, destroy, set

interface create
    module procedure :: create_all
    module procedure :: create_reals
    module procedure :: create_triangle_reals
end interface create

interface destroy
    module procedure :: destroy_all
end interface destroy

interface set
    module procedure :: set_energies
    module procedure :: add_single
    module procedure :: add_double
end interface set

contains

    pure subroutine create_all(energies, num_components)
        type(Concrete_Energies), intent(out) :: energies
        integer, intent(in) :: num_components

        call create(energies%walls_energies, num_components)
        call create(energies%short_energies, num_components)
        call create(energies%field_energies, num_components)
        call create(energies%dipolar_energies, num_components)
    end subroutine create_all

    pure subroutine destroy_all(energies)
        type(Concrete_Energies), intent(inout) :: energies

        call destroy_triangle_reals(energies%dipolar_energies)
        call destroy_reals(energies%field_energies)
        call destroy_triangle_reals(energies%short_energies)
        call destroy_reals(energies%walls_energies)
    end subroutine destroy_all

    pure subroutine add_single(energies, deltas, i_actor)
        type(Concrete_Energies), intent(inout) :: energies
        type(Concrete_Single_Energies), intent(in) :: deltas
        integer, intent(in) :: i_actor

        energies%walls_energies(i_actor) = energies%walls_energies(i_actor) + deltas%walls
        call add_energies(energies%short_energies, deltas%short, i_actor)
        energies%field_energies(i_actor) = energies%field_energies(i_actor) + deltas%field
        call add_energies(energies%dipolar_energies, deltas%dipolar, i_actor)
        energies%dipolar_mixture_energy = energies%dipolar_mixture_energy + deltas%dipolar_mixture
    end subroutine add_single

    pure subroutine add_double(energies, deltas, ij_actors)
        type(Concrete_Energies), intent(inout) :: energies
        type(Concrete_Double_Energies), intent(in) :: deltas
        integer, intent(in) :: ij_actors(:)

        integer :: i

        do i = 1, size(ij_actors)
            energies%walls_energies(ij_actors(i)) = energies%walls_energies(ij_actors(i)) + deltas%&
                walls(i)
            call add_energies(energies%short_energies, deltas%short(:, i), ij_actors(i))
            energies%field_energies(ij_actors(i)) = energies%field_energies(ij_actors(i)) + deltas%&
                field(i)
            call add_energies(energies%dipolar_energies, deltas%dipolar(:, i),  ij_actors(i))
        end do
        energies%dipolar_mixture_energy = energies%dipolar_mixture_energy + deltas%dipolar_mixture
    end subroutine add_double

    pure subroutine add_energies(energies, deltas, i_actor)
        type(Reals_Line), intent(inout) :: energies(:)
        real(DP), intent(in) :: deltas(:)
        integer, intent(in) :: i_actor

        integer :: i_component, i_observable, j_observable

        do i_component = 1, size(deltas)
            j_observable = maxval([i_actor, i_component])
            i_observable = minval([i_actor, i_component])
            energies(j_observable)%line(i_observable) = energies(j_observable)%&
                line(i_observable) + deltas(i_component)
        end do
    end subroutine add_energies

    pure subroutine set_energies(target_energies, source_energies)
        type(Concrete_Energies), intent(inout) :: target_energies
        type(Concrete_Energies), intent(in) :: source_energies

        target_energies%walls_energies = source_energies%walls_energies
        call set_energies_triangle(target_energies%short_energies, source_energies%short_energies)
        target_energies%field_energies = source_energies%field_energies
        call set_energies_triangle(target_energies%dipolar_energies, source_energies%&
            dipolar_energies)
        target_energies%dipolar_mixture_energy = source_energies%dipolar_mixture_energy
    end subroutine set_energies

    pure subroutine set_energies_triangle(energies_target, energies_source)
        type(Reals_Line), intent(inout) :: energies_target(:)
        type(Reals_Line), intent(in) :: energies_source(:)

        integer :: i_component

        do i_component = 1, size(energies_target)
            energies_target(i_component)%line = energies_source(i_component)%line
        end do
    end subroutine set_energies_triangle

end module procedures_observables_energies_factory
