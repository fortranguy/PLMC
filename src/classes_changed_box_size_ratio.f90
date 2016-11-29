module classes_changed_box_size_ratio

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Changed_Box_Size_Ratio
    private
        real(DP) :: delta = 0._DP
    contains
        procedure :: set => Abstract_set
        procedure :: get_delta => Abstract_get_delta
        procedure(Abstract_get), deferred :: get
    end type Abstract_Changed_Box_Size_Ratio

    abstract interface

        !> Box size change such that:
        !> \[ \ln(V^\prime) = \ln(V) + \mathrm{rand}[-1/2, +1/2] \delta \]
        function Abstract_get(this) result(ratio)
        import :: DP, num_dimensions, Abstract_Changed_Box_Size_Ratio
            class(Abstract_Changed_Box_Size_Ratio), intent(in) :: this
            real(DP) :: ratio(num_dimensions)
        end function Abstract_get

    end interface

    type, extends(Abstract_Changed_Box_Size_Ratio), public :: XYZ_Changed_Box_Size_Ratio
    contains
        procedure :: get => XYZ_get
    end type XYZ_Changed_Box_Size_Ratio

    type, extends(Abstract_Changed_Box_Size_Ratio), public :: XY_Changed_Box_Size_Ratio
    contains
        procedure :: get => XY_get
    end type XY_Changed_Box_Size_Ratio

    type, extends(Abstract_Changed_Box_Size_Ratio), public :: Null_Changed_Box_Size_Ratio
    contains
        procedure :: set => Null_set
        procedure :: get_delta => Null_get_delta
        procedure :: get => Null_get
    end type Null_Changed_Box_Size_Ratio

contains

!implementation Abstract_Changed_Box_Size_Ratio

    subroutine Abstract_set(this, delta)
        class(Abstract_Changed_Box_Size_Ratio), intent(inout) :: this
        real(DP), intent(in) :: delta

        call check_positive("Abstract_Changed_Box_Size_Ratio: construct", "delta", delta)
        this%delta = delta
    end subroutine Abstract_set

    pure real(DP) function Abstract_get_delta(this) result(delta)
        class(Abstract_Changed_Box_Size_Ratio), intent(in) :: this

        delta = this%delta
    end function Abstract_get_delta

!end implementation Abstract_Changed_Box_Size_Ratio

!implementation XYZ_Changed_Box_Size_Ratio

    !> \[ \frac{V^\prime}{V} = e^{\mathrm{rand}[-1/2, +1/2] \delta} \]
    function XYZ_get(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(XYZ_Changed_Box_Size_Ratio), intent(in) :: this

        real(DP) :: rand

        call random_number(rand)
        ratio = exp((rand - 0.5_DP) * this%delta / 3._DP)
    end function XYZ_get

!end implementation XYZ_Changed_Box_Size_Ratio

!implementation XY_Changed_Box_Size_Ratio

    !> \[ \frac{S^\prime}{S} = e^{\mathrm{rand}[-1/2, +1/2] \delta} \]
    function XY_get(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(XY_Changed_Box_Size_Ratio), intent(in) :: this

        real(DP) :: rand

        call random_number(rand)
        ratio(1:2) = exp((rand - 0.5_DP) * this%delta / 2._DP)
        ratio(3) = 1._DP
    end function XY_get

!end implementation XY_Changed_Box_Size_Ratio

!implementation Null_Changed_Box_Size_Ratio

    subroutine Null_set(this, delta)
        class(Null_Changed_Box_Size_Ratio), intent(inout) :: this
        real(DP), intent(in) :: delta
    end subroutine Null_set

    pure real(DP) function Null_get_delta(this) result(delta)
        class(Null_Changed_Box_Size_Ratio), intent(in) :: this
        delta = 0._DP
    end function Null_get_delta

    function Null_get(this) result(ratio)
        real(DP) :: ratio(num_dimensions)
        class(Null_Changed_Box_Size_Ratio), intent(in) :: this
        ratio = 1._DP
    end function Null_get

!end implementation Null_Changed_Box_Size_Ratio

end module classes_changed_box_size_ratio
