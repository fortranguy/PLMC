module procedures_field_expression_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private
public :: plate_expression

contains

    !> \( \vec{b}(x, z) \), cf. [[class_field_expression:Plates_get]]:
    !> \[
    !>      b_x(x, z) = \ln \left[ \frac{(x + a/2)^2 + z^2}{(x - a/2)^2 + z^2} \right] \\
    !>      b_y(x, z) = 0 \\
    !>      b_z(x, z) = 2 \left[
    !>          \arctan \left( \frac{x + a/2}{z} \right) -
    !>          \arctan \left( \frac{x - a/2}{z} \right)
    !>      \right].
    !> \]
    pure function plate_expression(size_x, position_13)
        real(DP), dimension(num_dimensions) :: plate_expression
        real(DP), intent(in) :: size_x
        real(DP), dimension(2), intent(in) :: position_13

        plate_expression(1) = log(((position_13(1) + size_x/2._DP)**2 + position_13(2)**2) / &
                                  ((position_13(1) - size_x/2._DP)**2 + position_13(2)**2))
        plate_expression(2) = 0._DP
        plate_expression(3) = 2._DP * (atan((position_13(1) + size_x/2._DP) / position_13(2)) - &
                                       atan((position_13(1) - size_x/2._DP) / position_13(2)))
    end function plate_expression

end module procedures_field_expression_micro
