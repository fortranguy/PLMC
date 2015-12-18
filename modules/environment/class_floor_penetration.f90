module class_floor_penetration

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Floor_Penetration
    contains
        procedure(Abstract_meet), deferred :: meet
    end type Abstract_Floor_Penetration

    abstract interface

        pure subroutine Abstract_meet(this, overlap, shortest_vector_from_floor, &
            position_from_floor)
            import :: DP, num_dimensions, Abstract_Floor_Penetration
            class(Abstract_Floor_Penetration), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
            real(DP), intent(in) :: position_from_floor(num_dimensions)
        end subroutine Abstract_meet

    end interface

    type, extends(Abstract_Floor_Penetration), public :: Null_Floor_Penetration
    contains
        procedure :: meet => Null_meet
    end type Null_Floor_Penetration

    type, extends(Abstract_Floor_Penetration), public :: Flat_Floor_Penetration
    contains
        procedure :: meet => Flat_meet
    end type Flat_Floor_Penetration

contains

!implementation Null_Floor_Penetration

    pure subroutine Null_meet(this, overlap, shortest_vector_from_floor, position_from_floor)
        class(Null_Floor_Penetration), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: shortest_vector_from_floor(num_dimensions)
        real(DP), intent(in) :: position_from_floor(num_dimensions)
        shortest_vector_from_floor = 0._DP
        overlap = .false.
    end subroutine Null_meet

!end implementation Null_Floor_Penetration

!implementation Flat_Floor_Penetration

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

end module class_floor_penetration
