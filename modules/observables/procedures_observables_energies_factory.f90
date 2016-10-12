module procedures_observables_energies_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line
use types_observables_energies, only: Concrete_Single_Energies, Concrete_Double_Energies, &
    Concrete_Observables_Energies
use procedures_reals_factory, only: reals_create => create, reals_destroy => destroy

implicit none

private
public :: create, destroy, set

interface create
    module procedure :: create_line
    module procedure :: create_element
end interface create

interface destroy
    module procedure :: destroy_element
    module procedure :: destroy_line
end interface

interface set
    module procedure :: set_energies
    module procedure :: add_single
    module procedure :: add_double
end interface set

contains

    pure subroutine create_line(energies, num_boxes, num_components)
        type(Concrete_Observables_Energies), allocatable, intent(out) :: energies(:)
        integer, intent(in) :: num_boxes, num_components

        integer :: i_box

        allocate(energies(num_boxes))
        do i_box = 1, size(energies)
            call create(energies(i_box), num_components)
        end do
    end subroutine create_line

    pure subroutine destroy_line(energies)
        type(Concrete_Observables_Energies), allocatable, intent(inout) :: energies(:)

        integer :: i_box

        if (allocated(energies)) then
            do i_box = size(energies), 1, -1
                call destroy(energies(i_box))
            end do
            deallocate(energies)
        end if
    end subroutine destroy_line

    pure subroutine create_element(energies, num_components)
        type(Concrete_Observables_Energies), intent(out) :: energies
        integer, intent(in) :: num_components

        allocate(energies%walls_energies(num_components))
        energies%walls_energies = 0._DP
        call reals_create(energies%short_energies, num_components)
        allocate(energies%field_energies(num_components))
        energies%field_energies = 0._DP
        call reals_create(energies%dipolar_energies, num_components)
    end subroutine create_element

    pure subroutine destroy_element(energies)
        type(Concrete_Observables_Energies), intent(inout) :: energies

        call reals_destroy(energies%dipolar_energies)
        if (allocated(energies%field_energies)) deallocate(energies%field_energies)
        call reals_destroy(energies%short_energies)
        if (allocated(energies%walls_energies)) deallocate(energies%walls_energies)
    end subroutine destroy_element

    pure subroutine add_single(energies, deltas, i_actor)
        type(Concrete_Observables_Energies), intent(inout) :: energies
        type(Concrete_Single_Energies), intent(in) :: deltas
        integer, intent(in) :: i_actor

        energies%walls_energies(i_actor) = energies%walls_energies(i_actor) + deltas%walls_energy
        call add_energies(energies%short_energies, deltas%short_energies, i_actor)
        energies%field_energies(i_actor) = energies%field_energies(i_actor) + deltas%field_energy
        call add_energies(energies%dipolar_energies, deltas%dipolar_energies, i_actor)
        energies%dipolar_shared_energy = energies%dipolar_shared_energy + deltas%&
            dipolar_shared_energy
    end subroutine add_single

    pure subroutine add_double(energies, deltas, ij_actors)
        type(Concrete_Observables_Energies), intent(inout) :: energies
        type(Concrete_Double_Energies), intent(in) :: deltas
        integer, intent(in) :: ij_actors(:)

        integer :: i

        do i = 1, size(ij_actors)
            energies%walls_energies(ij_actors(i)) = energies%walls_energies(ij_actors(i)) + deltas%&
                walls_energies(i)
            call add_energies(energies%short_energies, deltas%short_energies(:, i), ij_actors(i))
            energies%field_energies(ij_actors(i)) = energies%field_energies(ij_actors(i)) + deltas%&
                field_energies(i)
            call add_energies(energies%dipolar_energies, deltas%dipolar_energies(:, i),  &
                ij_actors(i))
        end do
        energies%dipolar_shared_energy = energies%dipolar_shared_energy + deltas%&
            dipolar_shared_energy
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
        type(Concrete_Observables_Energies), intent(inout) :: target_energies
        type(Concrete_Observables_Energies), intent(in) :: source_energies

        target_energies%walls_energies = source_energies%walls_energies
        call set_energies_triangle(target_energies%short_energies, source_energies%short_energies)
        target_energies%field_energies = source_energies%field_energies
        call set_energies_triangle(target_energies%dipolar_energies, source_energies%&
            dipolar_energies)
        target_energies%dipolar_shared_energy = source_energies%dipolar_shared_energy
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
