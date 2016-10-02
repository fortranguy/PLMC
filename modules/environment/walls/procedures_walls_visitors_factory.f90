module procedures_walls_visitors_factory

use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_walls_visitor, only: Abstract_Walls_Visitor, Concrete_Walls_Visitor, Null_Walls_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitors, visitable_walls, interact)
        class(Abstract_Walls_Visitor), allocatable, intent(out) :: visitors(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls(:)
        logical, intent(in) :: interact

        integer :: i_box

        if (interact) then
            allocate(Concrete_Walls_Visitor :: visitors(size(visitable_walls)))
        else
            allocate(Null_Walls_Visitor :: visitors(size(visitable_walls)))
        end if

        do i_box = 1, size(visitors)
            call visitors(i_box)%construct(visitable_walls(i_box))
        end do
    end subroutine create

    subroutine destroy(visitors)
        class(Abstract_Walls_Visitor), allocatable, intent(inout) :: visitors(:)

        integer :: i_box

        if (allocated(visitors)) then
            do i_box = size(visitors), 1, -1
                call visitors(i_box)%destroy()
            end do
            deallocate(visitors)
        end if
    end subroutine destroy

end module procedures_walls_visitors_factory
