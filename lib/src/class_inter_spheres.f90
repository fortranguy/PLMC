module class_inter_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_sphere, only: Abstract_Sphere

implicit none

private

    type, abstract, public :: Abstract_Inter_Spheres
    private
        real(DP) :: diameter
    contains
        procedure :: set_diameter => Abstract_Inter_Spheres_set_diameter
        procedure :: get_diameter => Abstract_Inter_Spheres_get_diameter
    end type Abstract_Inter_Spheres

    type, extends(Abstract_Inter_Spheres), public :: Concrete_Inter_Spheres

    end type Concrete_Inter_Spheres

contains

    subroutine Abstract_Inter_Spheres_set_diameter(this, sphere1, sphere2, non_additivity)
        class(Abstract_Inter_Spheres), intent(out) :: this
        class(Abstract_Sphere), intent(in) :: sphere1, sphere2
        real(DP), intent(in) :: non_additivity

        real(DP) :: diameter

        diameter = (sphere1%get_diameter() + sphere2%get_diameter())/2._DP + non_additivity
        if (diameter < 0._DP) call error_exit("Diameter is negative.")
        if (diameter < real_zero) call warning_continue("Diameter may be too small.")
        this%diameter = diameter
    end subroutine Abstract_Inter_Spheres_set_diameter

    pure function Abstract_Inter_Spheres_get_diameter(this) result(get_diameter)
        class(Abstract_Inter_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Abstract_Inter_Spheres_get_diameter

end module class_inter_spheres