module classes_min_distance

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Min_Distance
    private
        real(DP) :: min_distance = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Min_Distance

    type, extends(Abstract_Min_Distance), public :: Concrete_Min_Distance

    end type Concrete_Min_Distance

    type, extends(Abstract_Min_Distance), public :: Null_Min_Distance
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Min_Distance

    type, public :: Min_Distance_Wrapper
        class(Abstract_Min_Distance), allocatable :: distance
    end type Min_Distance_Wrapper

    type, public :: Min_Distance_Line
        type(Min_Distance_Wrapper), allocatable :: line(:)
    end type Min_Distance_Line

contains

!implementation Abstract_Min_Distance

    subroutine Abstract_set(this, min_distance)
        class(Abstract_Min_Distance), intent(inout) :: this
        real(DP), intent(in) :: min_distance

        call check_positive("Abstract_Min_Distance: set", "min_distance", min_distance)
        this%min_distance = min_distance
    end subroutine Abstract_set

    pure real(DP) function Abstract_get(this) result(min_distance)
        class(Abstract_Min_Distance), intent(in) :: this

        min_distance = this%min_distance
    end function Abstract_get

!end implementation Abstract_Min_Distance

!implementation Null_Min_Distance

    subroutine Null_set(this, min_distance)
        class(Null_Min_Distance), intent(inout) :: this
        real(DP), intent(in) :: min_distance
    end subroutine Null_set

    pure real(DP) function Null_get(this) result(min_distance)
        class(Null_Min_Distance), intent(in) :: this
        min_distance = 0._DP
    end function Null_get

!end implementation Null_Min_Distance

end module classes_min_distance
