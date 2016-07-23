module classes_dirac_distribution_plus

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Dirac_Distribution_Plus
    contains
        procedure(Abstrac_get_width), deferred :: get_width
        procedure(Abstract_get), deferred :: get
    end type Abstract_Dirac_Distribution_Plus

    abstract interface

        pure real(DP) function Abstrac_get_width(this)
            import :: DP, Abstract_Dirac_Distribution_Plus
            class(Abstract_Dirac_Distribution_Plus), intent(in) :: this
        end function Abstrac_get_width

        pure real(DP) function Abstract_get(this, distance)
        import :: DP, Abstract_Dirac_Distribution_Plus
            class(Abstract_Dirac_Distribution_Plus), intent(in) :: this
            real(DP), intent(in) :: distance
        end function Abstract_get

    end interface

    type, extends(Abstract_Dirac_Distribution_Plus), public :: Rectangular_Dirac_Distribution_Plus
    private
        real(DP) :: delta_distance = 0._DP !volume dependency?
    contains
        procedure :: set => Rectangular_set
        procedure :: get_width => Rectangular_get_width
        procedure :: get => Rectangular_get
    end type Rectangular_Dirac_Distribution_Plus

    type, extends(Abstract_Dirac_Distribution_Plus), public :: Gaussian_Dirac_Distribution_Plus
    private
        real(DP) :: max_distance = 0._DP
        real(DP) :: std_dev = 0._DP !! \( \sigma \)
    contains
        procedure :: set => Gaussian_set
        procedure :: get_width => Gaussian_get_width
        procedure :: get => Gaussian_get
    end type Gaussian_Dirac_Distribution_Plus

    type, extends(Abstract_Dirac_Distribution_Plus), public :: Null_Dirac_Distribution_Plus
    contains
        procedure :: get_width => Null_get_width
        procedure :: get => Null_get
    end type Null_Dirac_Distribution_Plus

contains

!implementation Rectangular_Dirac_Distribution_Plus

    subroutine Rectangular_set(this, delta_distance)
        class(Rectangular_Dirac_Distribution_Plus), intent(inout) :: this
        real(DP), intent(in) :: delta_distance

        call check_positive("Rectangular_Dirac_Distribution_Plus: set", "delta_distance", &
            delta_distance)
        this%delta_distance = delta_distance
    end subroutine Rectangular_set

    !> \[ \mathrm{d}r \]
    pure real(DP) function Rectangular_get_width(this) result(width)
        class(Rectangular_Dirac_Distribution_Plus), intent(in) :: this

        width = this%delta_distance
    end function Rectangular_get_width

    !> \[ \frac{1}{\mathrm{d}r} \]
    pure real(DP) function Rectangular_get(this, distance) result(distribution) !change name?
        class(Rectangular_Dirac_Distribution_Plus), intent(in) :: this
        real(DP), intent(in) :: distance

        if (distance <= this%delta_distance) then
            distribution = 1._DP/this%delta_distance
        else
            distribution = 0._DP
        end if
    end function Rectangular_get

!end implementation Rectangular_Dirac_Distribution_Plus

!implementation Gaussian_Dirac_Distribution_Plus

    subroutine Gaussian_set(this, max_distance, num_std_devs)
        class(Gaussian_Dirac_Distribution_Plus), intent(inout) :: this
        real(DP), intent(in) :: max_distance
        integer :: num_std_devs !! \( \mathsf{n}_\sigma \)

        call check_positive("Gaussian_Dirac_Distribution_Plus: set", "max_distance", max_distance)
        this%max_distance = max_distance
        call check_positive("Gaussian_Dirac_Distribution_Plus: set", "num_std_devs", num_std_devs)
        this%std_dev = this%max_distance / real(num_std_devs, DP)
    end subroutine Gaussian_set

    !> \[ \mathsf{n}_\sigma \sigma \]
    pure real(DP) function Gaussian_get_width(this) result(width)
        class(Gaussian_Dirac_Distribution_Plus), intent(in) :: this

        width = this%max_distance
    end function Gaussian_get_width

    !> \[
    !>      \frac{1}{\sigma \sqrt{2\pi}}
    !>          \exp\left( -\frac{\left( r - \frac{\mathsf{n}_\sigma}{2}\sigma \right)^2}
    !>          {2 \sigma^2} \right)
    !> \]
    pure real(DP) function Gaussian_get(this, distance) result(distribution)
        class(Gaussian_Dirac_Distribution_Plus), intent(in) :: this
        real(DP), intent(in) :: distance

        distribution = 1._DP / this%std_dev/sqrt(2._DP*PI) * &
            exp(-(distance - this%max_distance/2._DP)**2 / 2._DP/this%std_dev**2)
    end function Gaussian_get

!end implementation Gaussian_Dirac_Distribution_Plus

!implementation Null_Dirac_Distribution_Plus

    pure real(DP) function Null_get_width(this) result(width)
        class(Null_Dirac_Distribution_Plus), intent(in) :: this
        width = 0._DP
    end function Null_get_width

    pure real(DP) function Null_get(this, distance) result(distribution)
        class(Null_Dirac_Distribution_Plus), intent(in) :: this
        real(DP), intent(in) :: distance
        distribution = 0._DP
    end function Null_get

!end implementation Null_Dirac_Distribution_Plus

end module classes_dirac_distribution_plus
