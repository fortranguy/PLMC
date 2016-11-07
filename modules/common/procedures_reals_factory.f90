module procedures_reals_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_real_wrapper, only: Real_Line, Real_Triangle_Line

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_square
    module procedure :: create_triangle
end interface create

interface destroy
    module procedure :: destroy_triangle
    module procedure :: destroy_square
end interface destroy

contains

    pure subroutine create_square(square, num_elements_1, num_elements_2)
        type(Real_Triangle_Line), allocatable, intent(out) :: square(:)
        integer, intent(in) :: num_elements_1, num_elements_2

        integer :: i_element, j_element

        allocate(square(num_elements_1))
        do j_element = 1, size(square)
            allocate(square(j_element)%line(j_element))
            do i_element = 1, size(square(j_element)%line)
                call create(square(j_element)%line(i_element)%triangle, num_elements_2)
            end do
        end do
    end subroutine create_square

    pure subroutine destroy_square(square)
        type(Real_Triangle_Line), allocatable, intent(inout) :: square(:)

        integer :: i_element, j_element

        if (allocated(square)) then
            do j_element = size(square), 1, -1
                if (allocated(square(j_element)%line)) then
                    do i_element = size(square(j_element)%line), 1, -1
                        call destroy(square(j_element)%line(i_element)%triangle)
                    end do
                    deallocate(square(j_element)%line)
                end if
            end do
            deallocate(square)
        end if
    end subroutine destroy_square

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
