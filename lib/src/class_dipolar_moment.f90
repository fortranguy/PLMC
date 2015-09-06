module class_dipolar_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Moment
    private
        real(DP) :: norm
    contains
        procedure :: set_norm => Abstract_Dipolar_Moment_set_norm
        procedure :: get_norm => Abstract_Dipolar_Moment_get_norm
    end type Abstract_Dipolar_Moment

    type, extends(Abstract_Dipolar_Moment), public :: Concrete_Dipolar_Moment

    end type Concrete_Dipolar_Moment

contains

    subroutine Abstract_Dipolar_Moment_set_norm(this, norm)
        class(Abstract_Dipolar_Moment), intent(out) :: this
        real(DP), intent(in) :: norm

        if (norm < 0._DP) call error_exit("Dipolar moment is negative.")
        if (norm < real_zero) call warning_continue("Dipolar moment may be too small.")
        this%norm = norm
    end subroutine Abstract_Dipolar_Moment_set_norm

    pure function Abstract_Dipolar_Moment_get_norm(this) result(get_norm)
        class(Abstract_Dipolar_Moment), intent(in) :: this
        real(DP) :: get_norm
        
        get_norm = this%norm
    end function Abstract_Dipolar_Moment_get_norm

end module class_dipolar_moment