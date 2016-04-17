module class_component_chemical_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Component_Chemical_Potential
    private
        real(DP) :: density = 0._DP ! right place?
        real(DP) :: excess = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: get_density => Abstract_get_density
        procedure :: get_excess => Abstract_get_excess
    end type Abstract_Component_Chemical_Potential

    type, extends(Abstract_Component_Chemical_Potential), public :: &
        Concrete_Component_Chemical_Potential

    end type Concrete_Component_Chemical_Potential

    type, extends(Abstract_Component_Chemical_Potential), public :: &
        Null_Component_Chemical_Potential
    contains
        procedure :: set => Null_set
        procedure :: get_density => Null_get_density
        procedure :: get_excess => Null_get_excess
    end type Null_Component_Chemical_Potential

contains

!implementation Abstract_Component_Chemical_Potential

    subroutine Abstract_set(this, density, excess)
        class(Abstract_Component_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, excess

        call check_positive("Abstract_Component_Chemical_Potential", "density", density)
        this%density = density
        call check_positive("Abstract_Component_Chemical_Potential", "excess", excess)
        this%excess = excess
    end subroutine Abstract_set

    pure function Abstract_get_density(this) result(density)
        class(Abstract_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: density

        density = this%density
    end function Abstract_get_density

    pure function Abstract_get_excess(this) result(excess)
        class(Abstract_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: excess

        excess = this%excess
    end function Abstract_get_excess

!end implementation Abstract_Component_Chemical_Potential

!implementation Null_Component_Chemical_Potential

    subroutine Null_set(this, density, excess)
        class(Null_Component_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, excess
    end subroutine Null_set

    pure function Null_get_density(this) result(density)
        class(Null_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: density
        density = 0._DP
    end function Null_get_density

    pure function Null_get_excess(this) result(excess)
        class(Null_Component_Chemical_Potential), intent(in) :: this
        real(DP) :: excess
        excess = 0._DP
    end function Null_get_excess

!end implementation Null_Component_Chemical_Potential

end module class_component_chemical_potential
