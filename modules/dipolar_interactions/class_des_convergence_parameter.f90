module class_des_convergence_parameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_DES_Convergence_Parameter
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: alpha_x_box = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_box_edge => Abstract_get_box_edge
        procedure :: get => Abstract_get
    end type Abstract_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: &
        Concrete_DES_Convergence_Parameter
    end type Concrete_DES_Convergence_Parameter

    type, extends(Abstract_DES_Convergence_Parameter), public :: Null_DES_Convergence_Parameter
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_box_edge => Null_get_box_edge
        procedure :: get => Null_get
    end type Null_DES_Convergence_Parameter

contains

!implementation Abstract_DES_Convergence_Parameter

    subroutine Abstract_construct(this, periodic_box, alpha_x_box)
        class(Abstract_DES_Convergence_Parameter), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: alpha_x_box

        this%periodic_box => periodic_box
        call check_positive("Abstract_DES_Convergence_Parameter: construct", "alpha_x_box", &
            alpha_x_box)
        this%alpha_x_box = alpha_x_box
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Convergence_Parameter), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure real(DP) function Abstract_get_box_edge(this) result(box_edge)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: this

        real(DP) :: box_size(num_dimensions)

        box_size = this%periodic_box%get_size()
        box_edge = box_size(1)
    end function Abstract_get_box_edge

    pure real(DP) function Abstract_get(this) result(alpha)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: this

        alpha = this%alpha_x_box / this%get_box_edge()
    end function Abstract_get

!end implementation Abstract_DES_Convergence_Parameter

!implementation Null_DES_Convergence_Parameter

    subroutine Null_construct(this, periodic_box, alpha_x_box)
        class(Null_DES_Convergence_Parameter), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: alpha_x_box
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Convergence_Parameter), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_get_box_edge(this) result(box_edge)
        class(Null_DES_Convergence_Parameter), intent(in) :: this
        box_edge = 0._DP
    end function Null_get_box_edge

    pure real(DP) function Null_get(this) result(alpha)
        class(Null_DES_Convergence_Parameter), intent(in) :: this
        alpha = 0._DP
    end function Null_get

!end implementation Null_DES_Convergence_Parameter

end module class_des_convergence_parameter
