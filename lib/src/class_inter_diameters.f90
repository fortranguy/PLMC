module class_inter_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive
use class_diameters, only: Abstract_Diameters

implicit none

private

    type, abstract, public :: Abstract_Inter_Diameters
    private
        real(DP) :: non_additivity
        class(Abstract_Diameters), pointer :: diameters_1, diameters_2
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
        
        this%diameters_1 => diameters_1
        this%diameters_2 => diameters_2
        call check_positive("Abstract_Inter_Diameters", "non_additivity", non_additivity)
        this%non_additivity = non_additivity
    end subroutine Abstract_Inter_Diameters_construct
    
    subroutine Abstract_Inter_Diameters_destroy(this)
        class(Abstract_Inter_Diameters), intent(inout) :: this
        
        this%diameters_1 => null()
        this%diameters_2 => null()
    end subroutine Abstract_Inter_Diameters_destroy

    pure function Abstract_Inter_Diameters_get(this, i_particle, j_particle) result(diameter)
        class(Abstract_Inter_Diameters), intent(in) :: this
        integer, intent(in) :: i_particle, j_particle
        real(DP) :: diameter
        
        diameter = (this%diameters_1%get(i_particle) + this%diameters_2%get(j_particle)) / 2._DP + &
                   this%non_additivity
    end function Abstract_Inter_Diameters_get

end module class_inter_diameters
