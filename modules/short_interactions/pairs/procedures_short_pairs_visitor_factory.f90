module procedures_short_pairs_visitor_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor, Concrete_Short_Pairs_Visitor, &
    Null_Short_Pairs_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitor, periodic_box, interact)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Short_Pairs_Visitor :: visitor)
        else
            allocate(Null_Short_Pairs_Visitor :: visitor)
        end if
        call visitor%construct(periodic_box)
    end subroutine create

    subroutine destroy(visitor)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy

end module procedures_short_pairs_visitor_factory
