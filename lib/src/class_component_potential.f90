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
