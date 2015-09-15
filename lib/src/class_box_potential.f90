module class_box_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use module_particles, only: Concrete_Particle
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, public :: Box_Potential_Facade
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Pair_Potential), pointer :: pair_potential
    contains
        procedure :: construct => Box_Potential_Facade_construct
        procedure :: destroy => Box_Potential_Facade_destroy
        procedure, private :: set_positions => Box_Potential_Facade_set_positions
        procedure, private :: set_pair_potential => Box_Potential_Facade_set_pair_potential
        generic :: set => set_positions, set_pair_potential
        procedure :: visit => Box_Potential_Facade_visit
    end type Box_Potential_Facade

contains

    subroutine Box_Potential_Facade_construct(this, periodic_box)
        class(Box_Potential_Facade), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Box_Potential_Facade_construct

    subroutine Box_Potential_Facade_destroy(this)
        class(Box_Potential_Facade), intent(inout) :: this

        this%pair_potential => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Box_Potential_Facade_destroy

    subroutine Box_Potential_Facade_set_positions(this, positions)
        class(Box_Potential_Facade), intent(inout) :: this
        class(Abstract_Positions), target, intent(in) :: positions

        this%positions => positions
    end subroutine Box_Potential_Facade_set_positions

    subroutine Box_Potential_Facade_set_pair_potential(this, pair_potential)
        class(Box_Potential_Facade), intent(inout) :: this
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%pair_potential => pair_potential
    end subroutine Box_Potential_Facade_set_pair_potential

    pure subroutine Box_Potential_Facade_visit(this, same_type, particle, overlap, energy)
        class(Box_Potential_Facade), intent(inout) :: this
        logical, intent(in) :: same_type
        type(Concrete_Particle), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        real(DP) :: energy_i, distance
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, this%positions%get_num()
            if (same_type .and. particle%i_particle == i_particle) cycle
            distance = this%periodic_box%distance(particle%position, &
                                                  this%positions%get(i_particle))
            call this%pair_potential%meet(distance, overlap, energy_i)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Box_Potential_Facade_visit

end module class_box_potential
