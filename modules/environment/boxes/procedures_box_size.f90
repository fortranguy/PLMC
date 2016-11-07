module procedures_box_size

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: max_distance, edge

contains

    !> \[ r_\text{max} = \frac{\sqrt{3}}{2} L \]
    !> @note A 3D box is assumed.
    pure real(DP) function max_distance(box_size)
        real(DP), intent(in) :: box_size(:)

        max_distance = sqrt(3._DP) / 2._DP * box_size(1)
    end function max_distance

    !> @note This function seems silly. How to slice a result of a function (which is an array)?
    pure real(DP) function edge(box_size)
        real(DP), intent(in) :: box_size(:)

        edge = box_size(1)
    end function edge

end module procedures_box_size
