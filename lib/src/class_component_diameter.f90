module class_component_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Component_Diameter
    private
        real(DP) :: diameter
        real(DP) :: min_diameter
    contains
        procedure :: set => Abstract_Component_Diameter_set
        procedure :: get => Abstract_Component_Diameter_get
        procedure :: get_min => Abstract_Component_Diameter_get_min
    end type Abstract_Component_Diameter

    type, extends(Abstract_Component_Diameter), public :: Concrete_Component_Diameter

    end type Concrete_Component_Diameter

    type, extends(Abstract_Component_Diameter), public :: Null_Component_Diameter
    contains
        procedure :: set => Null_Component_Diameter_set
        procedure :: get => Null_Component_Diameter_get
        procedure :: get_min => Null_Component_Diameter_get_min
    end type Null_Component_Diameter

contains

!implementation Abstract_Component_Diameter

    subroutine Abstract_Component_Diameter_set(this, diameter, min_factor)
        class(Abstract_Component_Diameter), intent(inout) :: this
        real(DP), intent(in) :: diameter, min_factor

        call check_positive("Abstract_Component_Diameter", "diameter", diameter)
        this%diameter = diameter
        this%min_diameter = min_factor * this%diameter
    end subroutine Abstract_Component_Diameter_set

    pure function Abstract_Component_Diameter_get(this) result(diameter)
        class(Abstract_Component_Diameter), intent(in) :: this
        real(DP) :: diameter

        diameter = this%diameter
    end function Abstract_Component_Diameter_get

    pure function Abstract_Component_Diameter_get_min(this) result(min_diameter)
        class(Abstract_Component_Diameter), intent(in) :: this
        real(DP) :: min_diameter

        min_diameter = this%min_diameter
    end function Abstract_Component_Diameter_get_min

!end implementation Abstract_Component_Diameter

!implementation Null_Component_Diameter

    subroutine Null_Component_Diameter_set(this, diameter, min_factor)
        class(Null_Component_Diameter), intent(inout) :: this
        real(DP), intent(in) :: diameter, min_factor
    end subroutine Null_Component_Diameter_set

    pure function Null_Component_Diameter_get(this) result(diameter)
        class(Null_Component_Diameter), intent(in) :: this
        real(DP) :: diameter
        diameter = 0._DP
    end function Null_Component_Diameter_get

    pure function Null_Component_Diameter_get_min(this) result(min_diameter)
        class(Null_Component_Diameter), intent(in) :: this
        real(DP) :: min_diameter
        min_diameter = 0._DP
    end function Null_Component_Diameter_get_min

!end implementation Null_Component_Diameter

end module class_component_diameter
