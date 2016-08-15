module procedures_triangle_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private
public :: triangle_observables_init, triangle_observables_add, operator(-), &
    triangle_observables_sum

interface triangle_observables_add
    module procedure :: add_line
    module procedure :: add_triangle
end interface triangle_observables_add

interface operator(-)
    module procedure :: reals_diff
end interface

contains

    elemental subroutine triangle_observables_init(reals)
        type(Reals_Line), intent(inout) :: reals

        reals%line = 0._DP
    end subroutine triangle_observables_init

    pure subroutine add_line(reals, reals_i)
        type(Reals_Line), intent(inout) :: reals(:)
        real(DP), intent(in) :: reals_i(:)

        integer :: i_component
        do i_component = 1, size(reals)
            reals(i_component)%line(i_component) = reals(i_component)%line(i_component) + &
                reals_i(i_component)
        end do
    end subroutine add_line

    elemental subroutine add_triangle(reals, reals_i)
        type(Reals_Line), intent(inout) :: reals
        type(Reals_Line), intent(in) :: reals_i

        reals%line = reals%line + reals_i%line
    end subroutine add_triangle

    elemental function reals_diff(reals_1, reals_2)
        type(Reals_Line) :: reals_diff
        type(Reals_Line), intent(in) :: reals_1, reals_2

        allocate(reals_diff%line(size(reals_1%line)))
        reals_diff%line = reals_1%line - reals_2%line
    end function reals_diff

    pure real(DP) function triangle_observables_sum(reals)
        type(Reals_Line), intent(in) :: reals(:)

        triangle_observables_sum = sum(sum_line(reals))
    end function triangle_observables_sum

    elemental real(DP) function sum_line(reals)
        type(Reals_Line), intent(in) :: reals

        sum_line = sum(reals%line)
    end function sum_line

end module procedures_triangle_observables
