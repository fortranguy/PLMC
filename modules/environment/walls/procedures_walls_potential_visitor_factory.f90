module procedures_walls_potential_visitor_factory

use classes_walls_potential, only: Abstract_Walls_Potential
use classes_walls_potential_visitor, only: Abstract_Walls_Potential_Visitor, &
    Concrete_Walls_Potential_Visitor, Null_Walls_Potential_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, walls_potential, interact)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Walls_Potential_Visitor :: visitor)
        else
            allocate(Null_Walls_Potential_Visitor :: visitor)
        end if
        call visitor%construct(walls_potential)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_walls_potential_visitor_factory
