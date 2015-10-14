module class_component_moment_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Component_Moment_Norm
    private
        real(DP) :: norm
    contains
        procedure :: set => Abstract_Component_Moment_Norm_set
        procedure :: get => Abstract_Component_Moment_Norm_get
    end type Abstract_Component_Moment_Norm

    type, extends(Abstract_Component_Moment_Norm), public :: Null_Component_Moment_Norm
    contains
        procedure :: set => Null_Component_Moment_Norm_set
        procedure :: get => Null_Component_Moment_Norm_get
    end type Null_Component_Moment_Norm

    type, extends(Abstract_Component_Moment_Norm), public :: Concrete_Component_Moment_Norm

    end type Concrete_Component_Moment_Norm

contains

!implementation Abstract_Component_Moment_Norm

    subroutine Abstract_Component_Moment_Norm_set(this, norm)
        class(Abstract_Component_Moment_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm

        call check_positive("Abstract_Component_Moment_Norm", "norm", norm)
        this%norm = norm
    end subroutine Abstract_Component_Moment_Norm_set

    pure function Abstract_Component_Moment_Norm_get(this) result(norm)
        class(Abstract_Component_Moment_Norm), intent(in) :: this
        real(DP) :: norm

        norm = this%norm
    end function Abstract_Component_Moment_Norm_get

!end implementation Abstract_Component_Moment_Norm

!implementation Null_Component_Moment_Norm

    subroutine Null_Component_Moment_Norm_set(this, norm)
        class(Null_Component_Moment_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm
    end subroutine Null_Component_Moment_Norm_set

    pure function Null_Component_Moment_Norm_get(this) result(norm)
        class(Null_Component_Moment_Norm), intent(in) :: this
        real(DP) :: norm
        norm = 0._DP
    end function Null_Component_Moment_Norm_get

!end implementation Null_Component_Moment_Norm

end module class_component_moment_norm
