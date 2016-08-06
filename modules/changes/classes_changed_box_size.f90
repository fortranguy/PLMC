module classes_changed_box_size

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Changed_Box_Size
    private
        real(DP) :: delta
    contains
        procedure :: set => Abstract_set
        procedure(Abstract_get_ratio), deferred :: get_ratio
    end type Abstract_Changed_Box_Size

    abstract interface

        !> Box size change such that:
        !> \[ \frac{V^\prime}{V} = e^{\mathrm{rand}[-1/2, +1/2] \delta} \]
        function Abstract_get_ratio(this) result(ratio)
        import :: DP, num_dimensions, Abstract_Changed_Box_Size
            class(Abstract_Changed_Box_Size), intent(in) :: this
            real(DP) :: ratio(num_dimensions)
        end function Abstract_get_ratio

    end interface

    type, extends(Abstract_Changed_Box_Size), public :: XYZ_Changed_Box_Size
    contains
        procedure :: get_ratio => XYZ_get_ratio
    end type XYZ_Changed_Box_Size

    type, extends(Abstract_Changed_Box_Size), public :: XY_Changed_Box_Size
    contains
        procedure :: get_ratio => XY_get_ratio
    end type

contains

!implementation Abstract_Changed_Box_Size

    subroutine Abstract_set(this, delta)
        class(Abstract_Changed_Box_Size), intent(inout) :: this
        real(DP), intent(in) :: delta

        call check_positive("Abstract_Changed_Box_Size: set", "delta", delta)
        this%delta = delta
    end subroutine Abstract_set

!end implementation Abstract_Changed_Box_Size

!implementation XYZ_Changed_Box_Size

    !> \[ \frac{L_{1:3}^\prime}{L_{1:3}} = e^{\mathrm{rand}[-1/2, +1/2] \delta / 3} \]
    function XYZ_get_ratio(this) result(ratio)
        class(XYZ_Changed_Box_Size), intent(in) :: this
        real(DP) :: ratio(num_dimensions)

        real(DP) :: rand

        call random_number(rand)
        ratio = exp((rand - 0.5_DP) * this%delta / 3._DP)
    end function XYZ_get_ratio

!end implementation XYZ_Changed_Box_Size

!implementation XY_Changed_Box_Size

    !> \[
    !>      \left( \frac{L_{1:2}^\prime}{L_{1:2}}, \frac{L_3^\prime}{L_3} \right) =
    !>          \left( e^{\mathrm{rand}[-1/2, +1/2] \delta / 2}, 1 \right)
    !> \]
    function XY_get_ratio(this) result(ratio)
        class(XY_Changed_Box_Size), intent(in) :: this
        real(DP) :: ratio(num_dimensions)

        real(DP) :: rand

        call random_number(rand)
        ratio(1:2) = exp((rand - 0.5_DP) * this%delta / 2._DP)
        ratio(3) = 1._DP
    end function XY_get_ratio

!end implementation XY_Changed_Box_Size

end module classes_changed_box_size
