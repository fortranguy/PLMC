module class_particles_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Particles_Diameter
    private
        real(DP) :: diameter
        real(DP) :: min_diameter
    contains
        procedure :: set => Abstract_Particles_Diameter_set
        procedure :: get => Abstract_Particles_Diameter_get
        procedure :: get_min => Abstract_Particles_Diameter_get_min
    end type Abstract_Particles_Diameter

    type, extends(Abstract_Particles_Diameter), public :: Concrete_Particles_Diameter

    end type Concrete_Particles_Diameter

    type, extends(Abstract_Particles_Diameter), public :: Null_Particles_Diameter
    contains
        procedure :: set => Null_Particles_Diameter_set
        procedure :: get => Null_Particles_Diameter_get
        procedure :: get_min => Null_Particles_Diameter_get_min
    end type Null_Particles_Diameter

contains

!implementation Abstract_Particles_Diameter

    subroutine Abstract_Particles_Diameter_set(this, diameter, min_factor)
        class(Abstract_Particles_Diameter), intent(inout) :: this
        real(DP), intent(in) :: diameter, min_factor

        call check_positive("Abstract_Particles_Diameter", "diameter", diameter)
        this%diameter = diameter
        this%min_diameter = min_factor * this%diameter
    end subroutine Abstract_Particles_Diameter_set

    pure function Abstract_Particles_Diameter_get(this) result(diameter)
        class(Abstract_Particles_Diameter), intent(in) :: this
        real(DP) :: diameter

        diameter = this%diameter
    end function Abstract_Particles_Diameter_get

    pure function Abstract_Particles_Diameter_get_min(this) result(min_diameter)
        class(Abstract_Particles_Diameter), intent(in) :: this
        real(DP) :: min_diameter

        min_diameter = this%min_diameter
    end function Abstract_Particles_Diameter_get_min

!end implementation Abstract_Particles_Diameter

!implementation Null_Particles_Diameter

    subroutine Null_Particles_Diameter_set(this, diameter, min_factor)
        class(Null_Particles_Diameter), intent(inout) :: this
        real(DP), intent(in) :: diameter, min_factor
    end subroutine Null_Particles_Diameter_set

    pure function Null_Particles_Diameter_get(this) result(diameter)
        class(Null_Particles_Diameter), intent(in) :: this
        real(DP) :: diameter
        diameter = 0._DP
    end function Null_Particles_Diameter_get

    pure function Null_Particles_Diameter_get_min(this) result(min_diameter)
        class(Null_Particles_Diameter), intent(in) :: this
        real(DP) :: min_diameter
        min_diameter = 0._DP
    end function Null_Particles_Diameter_get_min

!end implementation Null_Particles_Diameter

end module class_particles_diameter
