module class_component_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Component_Potential
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: component_positions => null()
    contains
        procedure :: construct => Abstract_Component_Potential_construct
        procedure :: destroy => Abstract_Component_Potential_destroy
        procedure :: visit => Abstract_Component_Potential_visit
    end type Abstract_Component_Potential

    type, extends(Abstract_Component_Potential), public :: Concrete_Component_Potential

    end type Concrete_Component_Potential

    type, extends(Abstract_Component_Potential), public :: Null_Component_Potential
    contains
        procedure :: construct => Null_Component_Potential_construct
        procedure :: destroy => Null_Component_Potential_destroy
        procedure :: visit => Null_Component_Potential_visit
    end type Null_Component_Potential

contains

!implementation Abstract_Component_Potential

    subroutine Abstract_Component_Potential_construct(this, periodic_box, component_positions)
        class(Abstract_Component_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions

        this%periodic_box => periodic_box
        this%component_positions => component_positions
    end subroutine Abstract_Component_Potential_construct

    subroutine Abstract_Component_Potential_destroy(this)
        class(Abstract_Component_Potential), intent(inout) :: this

        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Component_Potential_destroy

    pure subroutine Abstract_Component_Potential_visit(this, overlap, energy, particle, &
        pair_potential, same_type)
        class(Abstract_Component_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(in) :: same_type

        real(DP) :: energy_i, distance
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, this%component_positions%get_num()
            if (same_type .and. particle%i == i_particle) cycle
            distance = this%periodic_box%distance(particle%position, &
                this%component_positions%get(i_particle))
            call pair_potential%meet(overlap, energy_i, distance)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Abstract_Component_Potential_visit

    pure subroutine Abstract_Component_Potential_visit_inter(this, overlap, energy, &
        component_1_positions, component_2_positions, pair_potential)
        class(Abstract_Component_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: component_1_positions, &
            component_2_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_ij, distance_ij
        integer :: i_particle, j_particle

        overlap = .false.
        energy = 0._DP
        do j_particle = 1, component_2_positions%get_num()
            do i_particle = 1, component_1_positions%get_num()
                distance_ij = this%periodic_box%distance(component_1_positions%get(i_particle), &
                    component_2_positions%get(j_particle))
                call pair_potential%meet(overlap, energy_ij, distance_ij)
                if (overlap) return
                energy = energy + energy_ij
            end do
        end do
    end subroutine Abstract_Component_Potential_visit_inter

    pure subroutine Abstract_Component_Potential_visit_intra(this, overlap, energy, &
        component_positions, pair_potential)
        class(Abstract_Component_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_ij, distance_ij
        integer :: i_particle, j_particle

        overlap = .false.
        energy = 0._DP
        do j_particle = 1, component_positions%get_num()
            do i_particle = j_particle + 1, component_positions%get_num()
                distance_ij = this%periodic_box%distance(component_positions%get(i_particle), &
                    component_positions%get(j_particle))
                call pair_potential%meet(overlap, energy_ij, distance_ij)
                if (overlap) return
                energy = energy + energy_ij
            end do
        end do
    end subroutine Abstract_Component_Potential_visit_intra

!end implementation Abstract_Component_Potential

!implementation Null_Component_Potential

    subroutine Null_Component_Potential_construct(this, periodic_box, component_positions)
        class(Null_Component_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
    end subroutine Null_Component_Potential_construct

    subroutine Null_Component_Potential_destroy(this)
        class(Null_Component_Potential), intent(inout) :: this
    end subroutine Null_Component_Potential_destroy

    pure subroutine Null_Component_Potential_visit(this, overlap, energy, particle, &
        pair_potential, same_type)
        class(Null_Component_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(in) :: same_type
        overlap = .false.
        energy = 0._DP
    end subroutine Null_Component_Potential_visit

!end implementation Null_Component_Potential

end module class_component_potential
