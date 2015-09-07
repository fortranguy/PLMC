module class_inter_spheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_spheres, only: Abstract_Spheres, Abstract_Spheres_Pointer

implicit none

private

    type, abstract, public :: Abstract_Inter_Spheres
    private
        type(Abstract_Spheres_Pointer) :: spheres_1, spheres_2
        real(DP) :: non_additivity
    contains
        procedure :: construct => Abstract_Inter_Spheres_construct
        procedure :: destroy => Abstract_Inter_Spheres_destroy
        procedure :: get_diameter => Abstract_Inter_Spheres_get_diameter
    end type Abstract_Inter_Spheres

    type, extends(Abstract_Inter_Spheres), public :: Concrete_Inter_Spheres

    end type Concrete_Inter_Spheres

contains

    subroutine Abstract_Inter_Spheres_construct(this, spheres_1, spheres_2, non_additivity)
        class(Abstract_Inter_Spheres), intent(out) :: this
        class(Abstract_Spheres), target, intent(in) :: spheres_1, spheres_2
        real(DP), intent(in) :: non_additivity
        
        this%spheres_1%ptr => spheres_1
        this%spheres_2%ptr => spheres_2
        if (non_additivity < 0._DP) then
            call error_exit("Abstract_Inter_Spheres: non_additivity is negative.")
        end if
        if (non_additivity < real_zero) then
            call warning_continue("Abstract_Inter_Spheres: non_additivity may be too small.")
        end if
        this%non_additivity = non_additivity
    end subroutine Abstract_Inter_Spheres_construct
    
    subroutine Abstract_Inter_Spheres_destroy(this)
        class(Abstract_Inter_Spheres), intent(inout) :: this
        
        this%spheres_1%ptr => null()
        this%spheres_2%ptr => null()
    end subroutine Abstract_Inter_Spheres_destroy

    pure function Abstract_Inter_Spheres_get_diameter(this, i_particle, j_particle) result(get_diameter)
        class(Abstract_Inter_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle, j_particle
        real(DP) :: get_diameter
        
        get_diameter = (this%spheres_1%ptr%get_diameter(i_particle) + &
                        this%spheres_2%ptr%get_diameter(j_particle)) / 2._DP + this%non_additivity
    end function Abstract_Inter_Spheres_get_diameter

end module class_inter_spheres
