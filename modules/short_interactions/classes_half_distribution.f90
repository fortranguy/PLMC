module classes_half_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Half_Distribution
    contains
        procedure(Abstrac_get_max_distance), deferred :: get_max_distance
        procedure(Abstract_get), deferred :: get
    end type Abstract_Half_Distribution

    abstract interface

        pure real(DP) function Abstrac_get_max_distance(this)
            import :: DP, Abstract_Half_Distribution
            class(Abstract_Half_Distribution), intent(in) :: this
        end function Abstrac_get_max_distance

        pure real(DP) function Abstract_get(this, distance)
        import :: DP, Abstract_Half_Distribution
            class(Abstract_Half_Distribution), intent(in) :: this
            real(DP), intent(in) :: distance
        end function Abstract_get

    end interface

    type, extends(Abstract_Half_Distribution), public :: Rectangular_Half_Distribution
    private
        real(DP) :: delta_distance !volume dependency?
    contains
        procedure :: set => Rectangular_set
        procedure :: get_max_distance => Rectangular_get_max_distance
        procedure :: get => Rectangular_get
    end type Rectangular_Half_Distribution

    type, extends(Abstract_Half_Distribution), public :: Null_Half_Distribution
    contains
        procedure :: get_max_distance => Null_get_max_distance
        procedure :: get => Null_get
    end type Null_Half_Distribution

contains

!implementation Rectangular_Half_Distribution

    subroutine Rectangular_set(this, delta_distance)
        class(Rectangular_Half_Distribution), intent(inout) :: this
        real(DP), intent(in) :: delta_distance

        call check_positive("Rectangular_Half_Distribution: set", "delta_distance", delta_distance)
        this%delta_distance = delta_distance
    end subroutine Rectangular_set

    pure real(DP) function Rectangular_get_max_distance(this) result(max_distance)
        class(Rectangular_Half_Distribution), intent(in) :: this

        max_distance = this%delta_distance
    end function Rectangular_get_max_distance

    pure real(DP) function Rectangular_get(this, distance) result(distribution) !change name?
        class(Rectangular_Half_Distribution), intent(in) :: this
        real(DP), intent(in) :: distance

        if (distance <= this%delta_distance) then
            distribution = 1._DP/this%delta_distance/2._DP
        else
            distribution = 0._DP
        end if
    end function Rectangular_get

!end implementation Rectangular_Half_Distribution

!implementation Null_Half_Distribution

    pure real(DP) function Null_get_max_distance(this) result(max_distance)
        class(Null_Half_Distribution), intent(in) :: this
        max_distance = 0._DP
    end function Null_get_max_distance

    pure real(DP) function Null_get(this, distance) result(distribution)
        class(Null_Half_Distribution), intent(in) :: this
        real(DP), intent(in) :: distance
        distribution = 0._DP
    end function Null_get

!end implementation Null_Half_Distribution

end module classes_half_distribution
