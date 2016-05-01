module classes_changed_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
implicit none

private

    type, abstract, public :: Abstract_Changed_Coordinates
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_increase_delta), deferred :: increase_delta
        procedure(Abstract_decrease_delta), deferred :: decrease_delta
        procedure(Abstract_get), deferred :: get
    end type Abstract_Changed_Coordinates

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_increase_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_increase_delta

        subroutine Abstract_decrease_delta(this)
        import :: Abstract_Changed_Coordinates
            class(Abstract_Changed_Coordinates), intent(inout) :: this
        end subroutine Abstract_decrease_delta

        function Abstract_get(this, i_particle) result(changed)
            import :: DP, num_dimensions, Abstract_Changed_Coordinates
            real(DP) :: changed(num_dimensions)
            class(Abstract_Changed_Coordinates), intent(in) :: this
            integer, intent(in) :: i_particle
        end function Abstract_get

    end interface

    type, extends(Abstract_Changed_Coordinates), public :: Null_Changed_Coordinates
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: increase_delta => Null_increase_delta
        procedure :: decrease_delta => Null_decrease_delta
        procedure :: get => Null_get
    end type Null_Changed_Coordinates

contains

!implementation Null_Changed_Coordinates

    subroutine Null_construct(this)
        class(Null_Changed_Coordinates), intent(out) :: this
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Changed_Coordinates), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_increase_delta(this)
        class(Null_Changed_Coordinates), intent(inout) :: this
    end subroutine Null_increase_delta

    subroutine Null_decrease_delta(this)
        class(Null_Changed_Coordinates), intent(inout) :: this
    end subroutine Null_decrease_delta

    function Null_get(this, i_particle) result(changed)
        real(DP) :: changed(num_dimensions)
        class(Null_Changed_Coordinates), intent(in) :: this
        integer, intent(in) :: i_particle
        changed = 0._DP
    end function Null_get

!end implementation Null_Changed_Coordinates

end module classes_changed_coordinates