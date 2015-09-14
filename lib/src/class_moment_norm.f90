module class_moment_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Moment_Norm
    private
        real(DP) :: norm
    contains
        procedure :: set => Abstract_Moment_Norm_set
        procedure :: get => Abstract_Moment_Norm_get
    end type Abstract_Moment_Norm

    type, extends(Abstract_Moment_Norm), public :: Null_Moment_Norm
    contains
        procedure :: set => Null_Moment_Norm_set
        procedure :: get => Null_Moment_Norm_get
    end type Null_Moment_Norm

    type, extends(Abstract_Moment_Norm), public :: Concrete_Moment_Norm

    end type Concrete_Moment_Norm

contains

!implementation Abstract_Moment_Norm

    subroutine Abstract_Moment_Norm_set(this, norm)
        class(Abstract_Moment_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm

        call check_positive("Abstract_Moment_Norm", "norm", norm)
        this%norm = norm
    end subroutine Abstract_Moment_Norm_set

    pure function Abstract_Moment_Norm_get(this) result(norm)
        class(Abstract_Moment_Norm), intent(in) :: this
        real(DP) :: norm

        norm = this%norm
    end function Abstract_Moment_Norm_get

!end implementation Abstract_Moment_Norm

!implementation Null_Moment_Norm

    subroutine Null_Moment_Norm_set(this, norm)
        class(Null_Moment_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm
    end subroutine Null_Moment_Norm_set

    pure function Null_Moment_Norm_get(this) result(norm)
        class(Null_Moment_Norm), intent(in) :: this
        real(DP) :: norm
        norm = 0._DP
    end function Null_Moment_Norm_get

!end implementation Null_Moment_Norm

end module class_moment_norm
