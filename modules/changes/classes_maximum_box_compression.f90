module classes_maximum_box_compression

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, abstract, public :: Abstract_Maximum_Box_Compression
    contains
        procedure(Abstract_get_delta), deferred, nopass :: get_delta
    end type Abstract_Maximum_Box_Compression

    abstract interface

        pure real(DP) function Abstract_get_delta(min_distance_ratio) result(delta)
        import :: DP
            real(DP), intent(in) :: min_distance_ratio
        end function Abstract_get_delta

    end interface

    type, extends(Abstract_Maximum_Box_Compression), public :: XYZ_Maximum_Box_Compression
    contains
        procedure, nopass :: get_delta => XYZ_get_delta
    end type XYZ_Maximum_Box_Compression

    type, extends(Abstract_Maximum_Box_Compression), public :: XY_Maximum_Box_Compression
    contains
        procedure, nopass :: get_delta => XY_get_delta
    end type XY_Maximum_Box_Compression

    type, extends(Abstract_Maximum_Box_Compression), public :: Null_Maximum_Box_Compression
    contains
        procedure, nopass :: get_delta => Null_get_delta
    end type Null_Maximum_Box_Compression

contains

!implementation XYZ_Maximum_Box_Compression

    !> \[ \delta = 6 \ln \left( \frac{r_\text{min}}{\sigma} \right) \]
    pure real(DP) function XYZ_get_delta(min_distance_ratio) result(delta)
        real(DP), intent(in) :: min_distance_ratio

        delta = 6._DP * log(min_distance_ratio)
    end function XYZ_get_delta

!end implementation XYZ_Maximum_Box_Compression

!implementation XY_Maximum_Box_Compression

    !> \[ \delta = 4 \ln \left( \frac{r_{1:2, \text{min}}}{\sigma} \right) \]
    pure real(DP) function XY_get_delta(min_distance_ratio) result(delta)
        real(DP), intent(in) :: min_distance_ratio

        delta = 4._DP * log(min_distance_ratio)
    end function XY_get_delta

!end implementation XY_Maximum_Box_Compression

!implementation Null_Maximum_Box_Compression

    pure real(DP) function Null_get_delta(min_distance_ratio) result(delta)
        real(DP), intent(in) :: min_distance_ratio
        delta = 0._DP
    end function Null_get_delta

!end implementation Null_Maximum_Box_Compression

end module classes_maximum_box_compression
