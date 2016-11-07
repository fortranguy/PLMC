module procedures_logical_factory

use types_logical_wrapper, only: Logical_Line

implicit none

private
public :: create

interface create
    module procedure :: create_triangle
end interface create

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

end module procedures_logical_factory
