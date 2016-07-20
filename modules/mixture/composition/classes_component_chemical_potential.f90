! remove from composition?
module classes_component_chemical_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Component_Chemical_Potential
    private
        real(DP) :: density = 0._DP ! right place?
        real(DP) :: inv_pow_activity = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: get_density => Abstract_get_density
        procedure :: get_inv_pow_activity => Abstract_get_inv_pow_activity
    end type Abstract_Component_Chemical_Potential

    type, extends(Abstract_Component_Chemical_Potential), public :: &
        Concrete_Component_Chemical_Potential

    end type Concrete_Component_Chemical_Potential

    type, extends(Abstract_Component_Chemical_Potential), public :: &
        Null_Component_Chemical_Potential
    contains
        procedure :: set => Null_set
        procedure :: get_density => Null_get_density
        procedure :: get_inv_pow_activity => Null_get_inv_pow_activity
    end type Null_Component_Chemical_Potential

contains

!implementation Abstract_Component_Chemical_Potential

    subroutine Abstract_set(this, density, inv_pow_activity)
        class(Abstract_Component_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, inv_pow_activity

        call check_positive("Abstract_Component_Chemical_Potential", "density", density)
        this%density = density
        call check_positive("Abstract_Component_Chemical_Potential", "inv_pow_activity", &
            inv_pow_activity)
        this%inv_pow_activity = inv_pow_activity
    end subroutine Abstract_set

    pure function Abstract_get_density(this) result(density)
        class(Abstract_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: density

        density = this%density
    end function Abstract_get_density

    !> \( a^{-N} \)
    pure function Abstract_get_inv_pow_activity(this) result(inv_pow_activity)
        class(Abstract_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: inv_pow_activity

        inv_pow_activity = this%inv_pow_activity
    end function Abstract_get_inv_pow_activity

!end implementation Abstract_Component_Chemical_Potential

!implementation Null_Component_Chemical_Potential

    subroutine Null_set(this, density, inv_pow_activity)
        class(Null_Component_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, inv_pow_activity
    end subroutine Null_set

    pure function Null_get_density(this) result(density)
        class(Null_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: density
        density = 0._DP
    end function Null_get_density

    pure function Null_get_inv_pow_activity(this) result(inv_pow_activity)
        class(Null_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: inv_pow_activity
        inv_pow_activity = 0._DP
    end function Null_get_inv_pow_activity

!end implementation Null_Component_Chemical_Potential

end module classes_component_chemical_potential
