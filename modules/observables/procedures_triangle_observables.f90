module procedures_triangle_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private
public :: triangle_observables_init, triangle_observables_add, triangle_observables_diff, &
    triangle_observables_sum

interface triangle_observables_add
    module procedure :: increment_line
    module procedure :: increment_triangle
end interface triangle_observables_add

interface triangle_observables_diff
    module procedure :: diff_triangle
end interface triangle_observables_diff

contains

    elemental subroutine triangle_observables_init(reals)
        type(Reals_Line), intent(inout) :: reals

        reals%line = 0._DP
    end subroutine triangle_observables_init

    pure subroutine increment_line(reals, reals_i)
        type(Reals_Line), intent(inout) :: reals(:)
        real(DP), intent(in) :: reals_i(:)

        integer :: i_component
        do i_component = 1, size(reals)
            reals(i_component)%line(i_component) = reals(i_component)%line(i_component) + &
                reals_i(i_component)
        end do
    end subroutine increment_line

    elemental subroutine increment_triangle(reals, reals_i)
        type(Reals_Line), intent(inout) :: reals
        type(Reals_Line), intent(in) :: reals_i

        reals%line = reals%line + reals_i%line
    end subroutine increment_triangle

    elemental subroutine diff_triangle(reals, reals_left, reals_right)
        type(Reals_Line), intent(inout) :: reals
        type(Reals_Line), intent(in) :: reals_left, reals_right

        reals%line = reals_left%line - reals_right%line
    end subroutine diff_triangle

    pure real(DP) function triangle_observables_sum(reals)
        type(Reals_Line), intent(in) :: reals(:)

        triangle_observables_sum = sum(sum_line(reals))
    end function triangle_observables_sum

    elemental real(DP) function sum_line(reals)
        type(Reals_Line), intent(in) :: reals

        sum_line = sum(reals%line)
    end function sum_line

end module procedures_triangle_observables
