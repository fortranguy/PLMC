module procedures_logical_factory

use types_logical_wrapper, only: Logical_Line, Logical_Rectangle

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_triangle
end interface create

interface destroy
    module procedure :: destroy_rectangles
end interface destroy

contains

    subroutine create_triangle(triangle, num_components)
        type(Logical_Line), allocatable, intent(out) :: triangle(:)
        integer, intent(in) :: num_components

        integer :: i_component

        allocate(triangle(num_components))
        do i_component = 1, size(triangle)
            allocate(triangle(i_component)%line(i_component))
            triangle(i_component)%line = .false.
        end do
    end subroutine create_triangle

    pure subroutine destroy_rectangles(rectangles)
        type(Logical_Rectangle), allocatable, intent(inout) ::rectangles(:, :, :)

        integer :: i_element, j_element, k_element

        if (allocated(rectangles)) then
            do k_element = size(rectangles, 3), 1, -1
                do j_element = size(rectangles, 2), 1, -1
                    do i_element = size(rectangles , 1), 1, -1
                        if (allocated(rectangles(i_element, j_element, k_element)%rectangle)) &
                            deallocate(rectangles(i_element, j_element, k_element)%rectangle)
                    end do
                end do
            end do
            deallocate(rectangles)
        end if
    end subroutine destroy_rectangles

end module procedures_logical_factory
