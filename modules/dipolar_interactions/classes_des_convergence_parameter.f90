module classes_des_convergence_parameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_DES_Convergence_Parameter
    private
        real(DP) :: alpha_x_box_edge = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: get_times_box_edge => Abstract_get_times_box_edge
    end type Abstract_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: &
        Concrete_DES_Convergence_Parameter
    end type Concrete_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: Null_DES_Convergence_Parameter
    contains
        procedure :: set => Null_set
        procedure :: get_times_box_edge => Null_get_times_box_edge
    end type Null_DES_Convergence_Parameter

contains

!implementation Abstract_DES_Convergence_Parameter

    subroutine Abstract_set(this, alpha_x_box_edge)
        class(Abstract_DES_Convergence_Parameter), intent(out) :: this
        real(DP), intent(in) :: alpha_x_box_edge

        call check_positive("Abstract_DES_Convergence_Parameter: set", "alpha_x_box_edge", &
            alpha_x_box_edge)
        this%alpha_x_box_edge = alpha_x_box_edge
    end subroutine Abstract_set

    pure real(DP) function Abstract_get_times_box_edge(this) result(alpha_x_box_edge)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: this

        alpha_x_box_edge = this%alpha_x_box_edge
    end function Abstract_get_times_box_edge

!end implementation Abstract_DES_Convergence_Parameter

!implementation Null_DES_Convergence_Parameter

    subroutine Null_set(this, alpha_x_box_edge)
        class(Null_DES_Convergence_Parameter), intent(out) :: this
        real(DP), intent(in) :: alpha_x_box_edge
    end subroutine Null_set

    pure real(DP) function Null_get_times_box_edge(this) result(alpha_x_box_edge)
        class(Null_DES_Convergence_Parameter), intent(in) :: this
        alpha_x_box_edge = 0._DP
    end function Null_get_times_box_edge

!end implementation Null_DES_Convergence_Parameter

end module classes_des_convergence_parameter
