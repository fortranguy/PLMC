module class_walls_potential_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_walls_potential, only: Abstract_Walls_Potential
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Walls_Potential_Visitor
    private
        class(Abstract_Walls_Potential), pointer :: walls_potential => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
    end type Abstract_Walls_Potential_Visitor

    type, extends(Abstract_Walls_Potential_Visitor), public :: Concrete_Walls_Potential_Visitor

    end type Concrete_Walls_Potential_Visitor

    type, extends(Abstract_Walls_Potential_Visitor), public :: Null_Walls_Potential_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
    end type Null_Walls_Potential_Visitor

contains

!implementation Abstract_Walls_Potential_Visitor

    subroutine Abstract_construct(this, walls_potential)
        class(Abstract_Walls_Potential_Visitor), intent(out) :: this
        class(Abstract_Walls_Potential), target, intent(in) :: walls_potential

        this%walls_potential => walls_potential
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Walls_Potential_Visitor), intent(inout) :: this

        this%walls_potential => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_visit(this, overlap, energy, component_positions, pair_potential)
        class(Abstract_Walls_Potential_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, component_positions%get_num()
            call this%walls_potential%visit(overlap, energy_i, component_positions%&
                get(i_particle), pair_potential)
            energy = energy + energy_i
        end do
    end subroutine Abstract_visit

!end implementation Abstract_Walls_Potential_Visitor

!implementation Null_Walls_Potential_Visitor

    subroutine Null_construct(this, walls_potential)
        class(Null_Walls_Potential_Visitor), intent(out) :: this
        class(Abstract_Walls_Potential), target, intent(in) :: walls_potential
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Walls_Potential_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit(this, overlap, energy, component_positions, pair_potential)
        class(Null_Walls_Potential_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit

!end implementation Null_Walls_Potential_Visitor

end module class_walls_potential_visitor
