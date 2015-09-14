module class_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Diameter
        real(DP) :: diameter
        real(DP) :: min_diameter
    contains
        procedure :: set => Abstract_Diameter_set
        procedure :: get => Abstract_Diameter_get
        procedure :: get_min => Abstract_Diameter_get_min
    end type Abstract_Diameter

    type, extends(Abstract_Diameter), public :: Concrete_Diameter

    end type Concrete_Diameter

contains

!implementation Abstract_Diameter

    subroutine Abstract_Diameter_set(this, diameter, min_factor)
        class(Abstract_Diameter), intent(inout) :: this
        real(DP), intent(in) :: diameter, min_factor

        call check_positive("Abstract_Diameter", "diameter", diameter)
        this%diameter = diameter
        this%min_diameter = min_factor * this%diameter
    end subroutine Abstract_Diameter_set

    pure function Abstract_Diameter_get(this) result(diameter)
        class(Abstract_Diameter), intent(in) :: this
        real(DP) :: diameter

        diameter = this%diameter
    end function Abstract_Diameter_get

    pure function Abstract_Diameter_get_min(this) result(min_diameter)
        class(Abstract_Diameter), intent(in) :: this
        real(DP) :: min_diameter

        min_diameter = this%min_diameter
    end function Abstract_Diameter_get_min

!end implementation Abstract_Diameter

end module class_diameter
