module procedures_centered_block_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private
public :: set_from_corner, set_from_wall

contains

    pure subroutine set_from_corner(shortest_vector, center, radius, position_13)
        real(DP), intent(out) :: shortest_vector(:)
        real(DP), intent(in) :: center(:), radius, position_13(:)

        real(DP) :: position_from_center(2)

        position_from_center = position_13 - center
        shortest_vector = (1._DP - radius/norm2(position_from_center)) * position_from_center
    end subroutine set_from_corner

    pure subroutine set_from_wall(shortest_vector, center, position_13)
        real(DP), intent(out) :: shortest_vector(:)
        real(DP), intent(in) :: center(:), position_13(:)

        real(DP) :: position_from_center(2)
        integer :: i_minloc

        position_from_center = position_13 - center
        i_minloc = minloc(position_from_center, 1)
        shortest_vector = 0._DP
        shortest_vector(i_minloc) = position_from_center(i_minloc)
    end subroutine set_from_wall

end module procedures_centered_block_micro
