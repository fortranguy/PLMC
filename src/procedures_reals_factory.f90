module procedures_reals_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_real_wrapper, only: Real_Line

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_triangle
end interface create

interface destroy
    module procedure :: destroy_triangle
end interface destroy

contains

    pure subroutine create_triangle(triangle, num_elements)
        type(Real_Line), allocatable, intent(out) :: triangle(:)
        integer, intent(in) :: num_elements

        integer :: i_element

        allocate(triangle(num_elements))
        do i_element = 1, size(triangle)
            allocate(triangle(i_element)%line(i_element))
            triangle(i_element)%line = 0._DP
        end do
    end subroutine create_triangle

    pure subroutine destroy_triangle(triangle)
        type(Real_Line), allocatable, intent(inout) :: triangle(:)

        integer :: i_element

        if (allocated(triangle)) then
            do i_element = size(triangle), 1, -1
                if (allocated(triangle(i_element)%line)) then
                    deallocate(triangle(i_element)%line)
                end if
            end do
            deallocate(triangle)
        end if
    end subroutine destroy_triangle

end module procedures_reals_factory
