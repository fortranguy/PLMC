module classes_moved_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
implicit none

private

    type, abstract, public :: Abstract_Moved_Coordinates
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_increase_delta), deferred :: increase_delta
        procedure(Abstract_decrease_delta), deferred :: decrease_delta
        procedure(Abstract_get), deferred :: get
    end type Abstract_Moved_Coordinates

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_increase_delta(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_increase_delta

        subroutine Abstract_decrease_delta(this)
        import :: Abstract_Moved_Coordinates
            class(Abstract_Moved_Coordinates), intent(inout) :: this
        end subroutine Abstract_decrease_delta

        function Abstract_get(this, i_particle) result(moved_coordinates)
            import :: DP, num_dimensions, Abstract_Moved_Coordinates
            real(DP) :: moved_coordinates(num_dimensions)
            class(Abstract_Moved_Coordinates), intent(in) :: this
            integer, intent(in) :: i_particle
        end function Abstract_get

    end interface

    type, extends(Abstract_Moved_Coordinates), public :: Null_Moved_Coordinates
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: increase_delta => Null_increase_delta
        procedure :: decrease_delta => Null_decrease_delta
        procedure :: get => Null_get
    end type Null_Moved_Coordinates

contains

!implementation Null_Moved_Coordinates

    subroutine Null_construct(this)
        class(Null_Moved_Coordinates), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Moved_Coordinates), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_increase_delta(this)
        class(Null_Moved_Coordinates), intent(inout) :: this
    end subroutine Null_increase_delta

    subroutine Null_decrease_delta(this)
        class(Null_Moved_Coordinates), intent(inout) :: this
    end subroutine Null_decrease_delta

    function Null_get(this, i_particle) result(moved_coordinates)
        real(DP) :: moved_coordinates(num_dimensions)
        class(Null_Moved_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        moved_coordinates = 0._DP
    end function Null_get

!end implementation Null_Moved_Coordinates

end module classes_moved_coordinates
