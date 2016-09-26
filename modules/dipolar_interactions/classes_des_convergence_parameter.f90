module classes_des_convergence_parameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use classes_box_size_memento, only: Abstract_Box_Size_Memento

implicit none

private

    type, abstract, public :: Abstract_DES_Convergence_Parameter
    private
        class(Abstract_Box_Size_Memento), pointer :: box_size_memento => null()
        real(DP) :: alpha_x_box_edge = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_times_box_edge => Abstract_get_times_box_edge
    end type Abstract_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: &
        Concrete_DES_Convergence_Parameter
    end type Concrete_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: Null_DES_Convergence_Parameter
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_times_box_edge => Null_get_times_box_edge
    end type Null_DES_Convergence_Parameter

contains

!implementation Abstract_DES_Convergence_Parameter

    subroutine Abstract_construct(this, alpha_x_box_edge)
        class(Abstract_DES_Convergence_Parameter), intent(out) :: this
        real(DP), intent(in) :: alpha_x_box_edge

        call check_positive("Abstract_DES_Convergence_Parameter: construct", "alpha_x_box_edge", &
            alpha_x_box_edge)
        this%alpha_x_box_edge = alpha_x_box_edge
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Convergence_Parameter), intent(inout) :: this

        this%box_size_memento => null()
    end subroutine Abstract_destroy

    pure real(DP) function Abstract_get_times_box_edge(this) result(alpha_x_box_edge)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: this

        alpha_x_box_edge = this%alpha_x_box_edge
    end function Abstract_get_times_box_edge

!end implementation Abstract_DES_Convergence_Parameter

!implementation Null_DES_Convergence_Parameter

    subroutine Null_construct(this, alpha_x_box_edge)
        class(Null_DES_Convergence_Parameter), intent(out) :: this
        real(DP), intent(in) :: alpha_x_box_edge
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Convergence_Parameter), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_get_times_box_edge(this) result(alpha_x_box_edge)
        class(Null_DES_Convergence_Parameter), intent(in) :: this
        alpha_x_box_edge = 0._DP
    end function Null_get_times_box_edge

!end implementation Null_DES_Convergence_Parameter

end module classes_des_convergence_parameter
