module class_minimum_distance

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Minimum_Distance
    private
        real(DP) :: min_distance
    contains
        procedure :: set => Abstract_set
        procedure :: get => Abstract_get
    end type Abstract_Minimum_Distance

    type, extends(Abstract_Minimum_Distance), public :: Concrete_Minimum_Distance

    end type Concrete_Minimum_Distance

    type, extends(Abstract_Minimum_Distance), public :: Null_Minimum_Distance
    contains
        procedure :: set => Null_set
        procedure :: get => Null_get
    end type Null_Minimum_Distance

contains

!implementation Abstract_Minimum_Distance

    subroutine Abstract_set(this, min_distance)
        class(Abstract_Minimum_Distance), intent(inout) :: this
        real(DP), intent(in) :: min_distance

        call check_positive("Abstract_Minimum_Distance", "min_distance", min_distance)
        this%min_distance = min_distance
    end subroutine Abstract_set

    pure function Abstract_get(this) result(min_distance)
        class(Abstract_Minimum_Distance), intent(in) :: this
        real(DP) :: min_distance

        min_distance = this%min_distance
    end function Abstract_get

!end implementation Abstract_Minimum_Distance

!implementation Null_Minimum_Distance

    subroutine Null_set(this, min_distance)
        class(Null_Minimum_Distance), intent(inout) :: this
        real(DP), intent(in) :: min_distance
    end subroutine Null_set

    pure function Null_get(this) result(min_distance)
        class(Null_Minimum_Distance), intent(in) :: this
        real(DP) :: min_distance
        min_distance = 0._DP
    end function Null_get

!end implementation Null_Minimum_Distance

end module class_minimum_distance
