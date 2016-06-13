module classes_floor_penetration

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_array_size, check_positive
use procedures_centered_block_micro, only: set_from_corner, set_from_wall

implicit none

private

    type, abstract, public :: Abstract_Floor_Penetration
    contains
        procedure(Abstract_get_min_depth), deferred :: get_min_depth
        procedure(Abstract_get_max_depth), deferred :: get_max_depth
        procedure(Abstract_meet), deferred :: meet
    end type Abstract_Floor_Penetration

    abstract interface

        pure real(DP) function Abstract_get_min_depth(this)
        import :: DP, Abstract_Floor_Penetration
            class(Abstract_Floor_Penetration), intent(in) :: this
        end function Abstract_get_min_depth

        pure real(DP) function Abstract_get_max_depth(this)
        import :: DP, Abstract_Floor_Penetration
            class(Abstract_Floor_Penetration), intent(in) :: this
        end function Abstract_get_max_depth

        pure subroutine Abstract_meet(this, overlap, shortest_vector_from_floor, &
            position_from_floor)
        import :: DP, num_dimensions, Abstract_Floor_Penetration
            class(Abstract_Floor_Penetration), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
            real(DP), intent(in) :: position_from_floor(num_dimensions)
        end subroutine Abstract_meet

    end interface

    type, extends(Abstract_Floor_Penetration), public :: Flat_Floor_Penetration
    contains
        procedure :: get_min_depth => Flat_get_depth
        procedure :: get_max_depth => Flat_get_depth
        procedure :: meet => Flat_meet
    end type Flat_Floor_Penetration

    !> This is a flat floor with a rounded block at the center, cf.
    !> modules/environment/centered_block_penetration.tex which shows the right half.
    !> When using [[Block_meet]], if a position is in a blue area,
    !> shortestVectorFromFloor's origin will be on a rounder corner. Otherwise (i.e. white area),
    !> it will be on a flat portion.
    type, extends(Abstract_Floor_Penetration), public :: Centered_Block_Penetration
    private
        real(DP), dimension(2) :: size = 0._DP
        real(DP) :: radius = 0._DP
        real(DP), dimension(2) :: lower_in = 0._DP, lower_out = 0._DP, upper_in = 0._DP, &
            upper_out = 0._DP! right centers
    contains
        procedure :: set => Block_set
        procedure :: get_min_depth => Block_get_min_depth
        procedure :: get_max_depth => Block_get_max_depth
        procedure :: meet => Block_meet
    end type Centered_Block_Penetration

    type, extends(Abstract_Floor_Penetration), public :: Null_Floor_Penetration
    contains
        procedure :: get_min_depth => Null_get_depth
        procedure :: get_max_depth => Null_get_depth
        procedure :: meet => Null_meet
    end type Null_Floor_Penetration

contains

!implementation Flat_Floor_Penetration

    pure real(DP) function Flat_get_depth(this) result(depth)
        class(Flat_Floor_Penetration), intent(in) :: this

        depth = 0._DP
    end function Flat_get_depth

    pure subroutine Flat_meet(this, overlap, shortest_vector_from_floor, position_from_floor)
        class(Flat_Floor_Penetration), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
        real(DP), intent(in) :: position_from_floor(num_dimensions)

        shortest_vector_from_floor = [0._DP, 0._DP, position_from_floor(3)]
        if (shortest_vector_from_floor(3) < 0._DP) then
            overlap = .true.
        else
            overlap = .false.
        end if
    end subroutine Flat_meet

!end implementation Flat_Floor_Penetration

!implementation Centered_Block_Penetration

    subroutine Block_set(this, size, radius)
        class(Centered_Block_Penetration), intent(out) :: this
        real(DP), intent(in) :: size(:), radius

        call check_array_size("Centered_Block_Penetration: set", "size", size, 2)
        call check_positive("Centered_Block_Penetration: set", "size", size)
        this%size = size
        call check_positive("Centered_Block_Penetration: set", "radius", radius)
        if (this%size(1) < 2.0_DP*radius) then
            call error_exit("Centered_Block_Penetration: set: 2*radius > size_x.")
        end if
        if (this%size(2) < 2.0_DP*radius) then
            call error_exit("Centered_Block_Penetration: set: 2*radius > size_z.")
        end if
        this%radius = radius

        this%upper_in = [this%size(1)/2._DP, this%size(2)]
        this%upper_out = [this%size(1)/2._DP - this%radius, this%size(2) - this%radius]
        this%lower_in  = [this%size(1)/2._DP + this%radius, this%radius]
        this%lower_out = [this%size(1)/2._DP, 0._DP]
    end subroutine Block_set

    pure real(DP) function Block_get_min_depth(this) result(min_depth)
        class(Centered_Block_Penetration), intent(in) :: this

        min_depth = this%size(2)
    end function Block_get_min_depth

    pure real(DP) function Block_get_max_depth(this) result(max_depth)
        class(Centered_Block_Penetration), intent(in) :: this

        max_depth = 0._DP
    end function Block_get_max_depth

    pure subroutine Block_meet(this, overlap, shortest_vector_from_floor, position_from_floor)
        class(Centered_Block_Penetration), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
        real(DP), intent(in) :: position_from_floor(num_dimensions)

        real(DP), dimension(2) :: shortest_vector, position_13

        if (0._DP < position_from_floor(1)) then
            position_13 = [+position_from_floor(1), position_from_floor(3)]
        else
            position_13 = [-position_from_floor(1), position_from_floor(3)]
        end if

        if (all(this%upper_out < position_13)) then !
            call set_from_corner(shortest_vector, this%upper_out, this%radius, position_13)
        else if (all(position_13 < this%lower_in)) then
            call set_from_corner(shortest_vector, this%lower_in, this%radius, position_13)
        else if (position_13(2) < this%size(2)/2._DP) then
            !> Frame: (0, \vec{e}_x, \vec{e}_z)
            call set_from_wall(shortest_vector, this%lower_out, position_13)
        else
            !> Frame: (0^\prime, -\vec{e}_x, -\vec{e}_z)
            call set_from_wall(shortest_vector, -this%upper_in, -position_13)
            shortest_vector = -shortest_vector
        end if

        if (any(shortest_vector < 0._DP )) then
            overlap = .true.
        else
            overlap = .false.
        end if
        shortest_vector_from_floor = [shortest_vector(1), 0._DP, shortest_vector(2)]
    end subroutine Block_meet

!implementation Centered_Block_Penetration

!implementation Null_Floor_Penetration

    pure real(DP) function Null_get_depth(this) result(depth)
        class(Null_Floor_Penetration), intent(in) :: this
        depth = 0._DP
    end function Null_get_depth

    pure subroutine Null_meet(this, overlap, shortest_vector_from_floor, position_from_floor)
        class(Null_Floor_Penetration), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
        real(DP), intent(in) :: position_from_floor(num_dimensions)
        shortest_vector_from_floor = 0._DP !Is it what I expect from a null object?
        overlap = .false.
    end subroutine Null_meet

!end implementation Null_Floor_Penetration

end module classes_floor_penetration
