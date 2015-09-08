module class_inter_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use module_error, only: warning_continue, error_exit
use class_diameters, only: Abstract_Diameters, Abstract_Diameters_Pointer

implicit none

private

    type, abstract, public :: Abstract_Inter_Diameters
    private
        type(Abstract_Diameters_Pointer) :: diameters_1, diameters_2
        real(DP) :: non_additivity
    contains
        procedure :: construct => Abstract_Inter_Diameters_construct
        procedure :: destroy => Abstract_Inter_Diameters_destroy
        procedure :: get => Abstract_Inter_Diameters_get
    end type Abstract_Inter_Diameters

    type, extends(Abstract_Inter_Diameters), public :: Concrete_Inter_Diameters

    end type Concrete_Inter_Diameters

contains

    subroutine Abstract_Inter_Diameters_construct(this, diameters_1, diameters_2, non_additivity)
        class(Abstract_Inter_Diameters), intent(out) :: this
        class(Abstract_Diameters), target, intent(in) :: diameters_1, diameters_2
        real(DP), intent(in) :: non_additivity
        
        this%diameters_1%ptr => diameters_1
        this%diameters_2%ptr => diameters_2
        if (non_additivity < 0._DP) then
            call error_exit("Abstract_Inter_Diameters: non_additivity is negative.")
        end if
        if (non_additivity < real_zero) then
            call warning_continue("Abstract_Inter_Diameters: non_additivity may be too small.")
        end if
        this%non_additivity = non_additivity
    end subroutine Abstract_Inter_Diameters_construct
    
    subroutine Abstract_Inter_Diameters_destroy(this)
        class(Abstract_Inter_Diameters), intent(inout) :: this
        
        this%diameters_1%ptr => null()
        this%diameters_2%ptr => null()
    end subroutine Abstract_Inter_Diameters_destroy

    pure function Abstract_Inter_Diameters_get(this, i_particle, j_particle) result(diameter)
        class(Abstract_Inter_Diameters), intent(in) :: this
        integer, intent(in) :: i_particle, j_particle
        real(DP) :: diameter
        
        diameter = (this%diameters_1%ptr%get(i_particle) + this%diameters_2%ptr%get(j_particle)) / &
                   2._DP + this%non_additivity
    end function Abstract_Inter_Diameters_get

end module class_inter_diameters
