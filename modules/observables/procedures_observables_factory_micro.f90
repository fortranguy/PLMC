module procedures_observables_factory_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private
public :: create_reals, create_triangle_reals, create_triangle_nodes, create_square_reals, &
    destroy_reals, destroy_triangle_reals, destroy_triangle_nodes, destroy_square_reals

contains

    pure subroutine create_reals(reals, num_components)
        real(DP), allocatable, intent(out) :: reals(:)
        integer, intent(in) :: num_components

        allocate(reals(num_components))
    end subroutine create_reals

    pure subroutine destroy_reals(reals)
        real(DP), allocatable, intent(inout) :: reals(:)

        if (allocated(reals)) deallocate(reals)
    end subroutine destroy_reals

    pure subroutine create_triangle_reals(triangle, num_components)
        type(Reals_Line), allocatable, intent(out) :: triangle(:)
        integer, intent(in) :: num_components

        allocate(triangle(num_components))
        call create_triangle_nodes(triangle)
    end subroutine create_triangle_reals

    pure subroutine create_square_reals(square, num_components)
        real(DP), allocatable, intent(out) :: square(:, :)
        integer, intent(in) :: num_components

        allocate(square(num_components, num_components))
        square = 0._DP
    end subroutine create_square_reals

    pure subroutine destroy_square_reals(square)
        real(DP), allocatable, intent(out) :: square(:, :)

        if (allocated(square)) deallocate(square)
    end subroutine destroy_square_reals

    pure subroutine destroy_triangle_reals(triangle)
        type(Reals_Line), allocatable, intent(inout) :: triangle(:)

        if (allocated(triangle)) then
            call destroy_triangle_nodes(triangle)
            deallocate(triangle)
        end if
    end subroutine destroy_triangle_reals

    pure subroutine create_triangle_nodes(triangle)
        type(Reals_Line), intent(out) :: triangle(:)

        integer :: i_component
        do i_component = 1, size(triangle)
            allocate(triangle(i_component)%line(i_component))
            triangle(i_component)%line = 0._DP
        end do
    end subroutine create_triangle_nodes

    pure subroutine destroy_triangle_nodes(triangle)
        type(Reals_Line), intent(inout) :: triangle(:)

        integer :: i_component

        do i_component = size(triangle), 1, -1
            if (allocated(triangle(i_component)%line)) then
                deallocate(triangle(i_component)%line)
            end if
        end do
    end subroutine destroy_triangle_nodes

end module procedures_observables_factory_micro
