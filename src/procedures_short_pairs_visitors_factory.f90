module procedures_short_pairs_visitors_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor, Concrete_Short_Pairs_Visitor, &
    Null_Short_Pairs_Visitor

implicit none

private
public :: create, destroy

contains

    subroutine create(visitors, periodic_boxes, interact)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(out) :: visitors(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        logical, intent(in) :: interact

        integer :: i_box

        if (interact) then
            allocate(Concrete_Short_Pairs_Visitor :: visitors(size(periodic_boxes)))
        else
            allocate(Null_Short_Pairs_Visitor :: visitors(size(periodic_boxes)))
        end if

        do i_box = 1, size(visitors)
            call visitors(i_box)%construct(periodic_boxes(i_box))
        end do
    end subroutine create

    subroutine destroy(visitors)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(inout) :: visitors(:)

        integer :: i_box

        if (allocated(visitors)) then
            do i_box = size(visitors), 1, -1
                call visitors(i_box)%destroy()
            end do
            deallocate(visitors)
        end if
    end subroutine destroy

end module procedures_short_pairs_visitors_factory
