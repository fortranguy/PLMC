module procedures_walls_visitor_factory

use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_walls_visitor, only: Abstract_Walls_Visitor, Concrete_Walls_Visitor, Null_Walls_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, visitable_walls, interact)
        class(Abstract_Walls_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Walls_Visitor :: visitor)
        else
            allocate(Null_Walls_Visitor :: visitor)
        end if
        call visitor%construct(visitable_walls)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_Walls_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_walls_visitor_factory
