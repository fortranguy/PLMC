module class_walls_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_floor_penetration, only: Abstract_Floor_Penetration
use class_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Walls_Potential
    private
        real(DP) :: gap
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: are_outside_box => Abstract_are_outside_box
        procedure :: get_gap => Abstract_get_gap
        procedure :: visit => Abstract_visit
        procedure, private :: position_from_floor => Abstract_position_from_floor
        procedure, private :: position_from_ceiling => Abstract_position_from_ceiling
    end type Abstract_Walls_Potential

    type, extends(Abstract_Walls_Potential), public :: Null_Walls_Potential
        contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: are_outside_box => Null_are_outside_box
        procedure :: get_gap => Null_get_gap
        procedure :: visit => Null_visit
    end type Null_Walls_Potential

    type, extends(Abstract_Walls_Potential), public :: Concrete_Walls_Potential

    end type Concrete_Walls_Potential

contains

!implementation Abstract_Walls_Potential

    subroutine Abstract_construct(this, periodic_box, gap, floor_penetration)
        class(Abstract_Walls_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration

        this%periodic_box => periodic_box
        call check_positive("Abstract_Walls_Potential", "gap", gap)
        this%gap = gap
        if (2._DP * floor_penetration%get_min_depth() > this%gap) then
            call error_exit("Abstract_Walls_Potential: floor_penetration's min_depth is too big.")
        end if
        allocate(this%floor_penetration, source = floor_penetration)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Walls_Potential), intent(inout) :: this

        if (allocated(this%floor_penetration)) deallocate(this%floor_penetration)
        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> This function tells if the walls are outside the box.
    !> They must remain inside if box_size(3) is constant (which is expected in XY periodicity).
    pure logical function Abstract_are_outside_box(this) result(are_outside)
        class(Abstract_Walls_Potential), intent(in) :: this

        real(DP) :: box_size(num_dimensions)

        box_size = this%periodic_box%get_size()
        if (this%gap > box_size(3)) then
            are_outside = .true.
        else
            are_outside = .false.
        end if
    end function Abstract_are_outside_box

    pure real(DP) function Abstract_get_gap(this) result(gap)
        class(Abstract_Walls_Potential), intent(in) :: this

        gap = this%gap
    end function Abstract_get_gap

    pure subroutine Abstract_visit(this, overlap, energy, position, pair_potential)
        class(Abstract_Walls_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_floor, energy_ceiling
        real(DP), dimension(num_dimensions) :: shortest_vector_from_floor, &
            shortest_vector_from_ceiling

        call this%floor_penetration%meet(overlap, shortest_vector_from_floor, &
            this%position_from_floor(position))
        if (overlap) return
        call this%floor_penetration%meet(overlap, shortest_vector_from_ceiling, &
            -this%position_from_ceiling(position))
        if (overlap) return
        call pair_potential%meet(overlap, energy_floor, norm2(shortest_vector_from_floor))
        if (overlap) return
        call pair_potential%meet(overlap, energy_ceiling, norm2(shortest_vector_from_ceiling))
        if (overlap) return
        energy = energy_floor + energy_ceiling
    end subroutine Abstract_visit

    pure function Abstract_position_from_floor(this, position) result(position_from_floor)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_floor(num_dimensions)

        position_from_floor = position + [0._DP, 0._DP, this%gap/2._DP]
    end function Abstract_position_from_floor

    pure function Abstract_position_from_ceiling(this, position) result(position_from_ceiling)
        class(Abstract_Walls_Potential), intent(in) :: this
        real(DP), intent(in) :: position(:)
        real(DP) :: position_from_ceiling(num_dimensions)

        position_from_ceiling = position - [0._DP, 0._DP, this%gap/2._DP]
    end function Abstract_position_from_ceiling

!end implementation Abstract_Walls_Potential

!implementation Null_Walls_Potential

    subroutine Null_construct(this, periodic_box, gap, floor_penetration)
        class(Null_Walls_Potential), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: gap
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
    end subroutine Null_construct

    pure logical function Null_are_outside_box(this) result(are_outside)
        class(Null_Walls_Potential), intent(in) :: this
        are_outside = .false.
    end function Null_are_outside_box

    pure real(DP) function Null_get_gap(this) result(gap)
        class(Null_Walls_Potential), intent(in) :: this
        gap = 0._DP
    end function Null_get_gap

    pure subroutine Null_visit(this, overlap, energy, position, pair_potential)
        class(Null_Walls_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: position(:)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit

    subroutine Null_destroy(this)
        class(Null_Walls_Potential), intent(inout) :: this
    end subroutine Null_destroy

!end implementation Null_Walls_Potential

end module class_walls_potential
