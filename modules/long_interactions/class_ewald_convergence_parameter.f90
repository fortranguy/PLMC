module class_ewald_convergence_parameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box

implicit none

private

    type, abstract, public :: Abstract_Ewald_Convergence_Parameter
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        real(DP) :: alpha_x_box
    contains
        procedure :: construct => Abstract_Ewald_Convergence_Parameter_construct
        procedure :: destroy => Abstract_Ewald_Convergence_Parameter_destroy
        procedure :: get => Abstract_Ewald_Convergence_Parameter_get
    end type Abstract_Ewald_Convergence_Parameter

    type, extends(Abstract_Ewald_Convergence_Parameter), public :: &
        Concrete_Ewald_Convergence_Parameter
    end type Concrete_Ewald_Convergence_Parameter

    type, extends(Abstract_Ewald_Convergence_Parameter), public :: Null_Ewald_Convergence_Parameter
    end type Null_Ewald_Convergence_Parameter

contains

!implementation Abstract_Ewald_Convergence_Parameter

    subroutine Abstract_Ewald_Convergence_Parameter_construct(this, periodic_box, alpha_x_box)
        class(Abstract_Ewald_Convergence_Parameter), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: alpha_x_box

        this%periodic_box => periodic_box
        call check_positive("Abstract_Ewald_Convergence_Parameter_construct", "alpha_x_box", &
            alpha_x_box)
        this%alpha_x_box = alpha_x_box
    end subroutine Abstract_Ewald_Convergence_Parameter_construct

    subroutine Abstract_Ewald_Convergence_Parameter_destroy(this)
        class(Abstract_Ewald_Convergence_Parameter), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_Ewald_Convergence_Parameter_destroy

    pure real(DP) function Abstract_Ewald_Convergence_Parameter_get(this) result(alpha)
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: this

        real(DP) :: box_size(num_dimensions)

        box_size = this%periodic_box%get_size()
        alpha = this%alpha_x_box / box_size(1)
    end function Abstract_Ewald_Convergence_Parameter_get

!end implementation Abstract_Ewald_Convergence_Parameter

!implementation Null_Ewald_Convergence_Parameter

    subroutine Null_Ewald_Convergence_Parameter_construct(this, periodic_box, alpha_x_box)
        class(Null_Ewald_Convergence_Parameter), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        real(DP), intent(in) :: alpha_x_box
    end subroutine Null_Ewald_Convergence_Parameter_construct

    subroutine Null_Ewald_Convergence_Parameter_destroy(this)
        class(Null_Ewald_Convergence_Parameter), intent(inout) :: this
    end subroutine Null_Ewald_Convergence_Parameter_destroy

    pure real(DP) function Null_Ewald_Convergence_Parameter_get(this) result(alpha)
        class(Null_Ewald_Convergence_Parameter), intent(in) :: this
        alpha = 0._DP
    end function Null_Ewald_Convergence_Parameter_get

!end implementation Null_Ewald_Convergence_Parameter

end module class_ewald_convergence_parameter
