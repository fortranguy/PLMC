module procedures_walls_potential_visitor_factory

use classes_walls_potential, only: Abstract_Walls_Potential
use classes_walls_potential_visitor, only: Abstract_Walls_Potential_Visitor, &
    Concrete_Walls_Potential_Visitor, Null_Walls_Potential_Visitor

implicit none

private
public :: walls_potential_visitor_create, walls_potential_visitor_destroy

contains

    subroutine walls_potential_visitor_create(visitor, walls_potential, interact)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Walls_Potential_Visitor :: visitor)
        else
            allocate(Null_Walls_Potential_Visitor :: visitor)
        end if
        call visitor%construct(walls_potential)
    end subroutine walls_potential_visitor_create

    subroutine walls_potential_visitor_destroy(visitor)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine walls_potential_visitor_destroy

end module procedures_walls_potential_visitor_factory
