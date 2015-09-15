module class_box_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_periodic_box, only: Abstract_Periodic_Box
use class_positions, only: Abstract_Positions
use module_particles, only: Concrete_Particle
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Box_Potential
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box
        class(Abstract_Positions), pointer :: positions
        class(Abstract_Pair_Potential), pointer :: pair_potential
    contains
        procedure :: construct => Abstract_Box_Potential_construct
        procedure :: destroy => Abstract_Box_Potential_destroy
        procedure :: visit => Abstract_Box_Potential_visit
    end type Abstract_Box_Potential

contains

    subroutine Abstract_Box_Potential_construct(this, periodic_box, positions, pair_potential)
        class(Abstract_Box_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Positions), target, intent(in) :: positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        this%periodic_box => periodic_box
        this%positions => positions
        this%pair_potential => pair_potential
    end subroutine Abstract_Box_Potential_construct

    subroutine Abstract_Box_Potential_destroy(this)
        class(Abstract_Box_Potential), intent(inout) :: this

        this%pair_potential => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Box_Potential_destroy

    pure subroutine Abstract_Box_Potential_visit(this, same_type, particle, overlap, energy)
        class(Abstract_Box_Potential), intent(inout) :: this
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
    end subroutine Abstract_Box_Potential_visit

end module class_box_potential
