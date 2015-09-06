module class_sphere

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit

implicit none

private

    type, abstract, public :: Abstract_Sphere
    private
        real(DP) :: diameter
    contains
        procedure :: set_diameter => Abstract_Sphere_set_diameter
        procedure :: get_diameter => Abstract_Sphere_get_diameter
    end type Abstract_Sphere

    type, extends(Abstract_Sphere), public :: Concrete_Sphere

    end type Concrete_Sphere

contains

    subroutine Abstract_Sphere_set_diameter(this, diameter)
        class(Abstract_Sphere), intent(out) :: this
        real(DP), intent(in) :: diameter

        if (diameter < 0._DP) call error_exit("Diameter is negative.")
        if (diameter < real_zero) call warning_continue("Diameter may be too small.")
        this%diameter = diameter
    end subroutine Abstract_Sphere_set_diameter

    pure function Abstract_Sphere_get_diameter(this) result(get_diameter)
        class(Abstract_Sphere), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Abstract_Sphere_get_diameter

end module class_sphere