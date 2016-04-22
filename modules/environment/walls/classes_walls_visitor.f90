module classes_walls_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Walls_Visitor
    private
        class(Abstract_Visitable_Walls), pointer :: walls => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
    end type Abstract_Walls_Visitor

    type, extends(Abstract_Walls_Visitor), public :: Concrete_Walls_Visitor

    end type Concrete_Walls_Visitor

    type, extends(Abstract_Walls_Visitor), public :: Null_Walls_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
    end type Null_Walls_Visitor

contains

!implementation Abstract_Walls_Visitor

    subroutine Abstract_construct(this, walls)
        class(Abstract_Walls_Visitor), intent(out) :: this
        class(Abstract_Visitable_Walls), target, intent(in) :: walls

        this%walls => walls
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Walls_Visitor), intent(inout) :: this

        this%walls => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_visit(this, overlap, energy, positions, pair_potential)
        class(Abstract_Walls_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_i
        integer :: i_particle

        overlap = .false.
        energy = 0._DP
        do i_particle = 1, positions%get_num()
            call this%walls%visit(overlap, energy_i, positions%get(i_particle), &
                pair_potential)
            if (overlap) return
            energy = energy + energy_i
        end do
    end subroutine Abstract_visit

!end implementation Abstract_Walls_Visitor

!implementation Null_Walls_Visitor

    subroutine Null_construct(this, walls)
        class(Null_Walls_Visitor), intent(out) :: this
        class(Abstract_Visitable_Walls), target, intent(in) :: walls
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Walls_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit(this, overlap, energy, positions, pair_potential)
        class(Null_Walls_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit

!end implementation Null_Walls_Visitor

end module classes_walls_visitor
