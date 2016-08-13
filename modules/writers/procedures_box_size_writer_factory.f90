module procedures_box_size_writer_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_writer, only: Abstract_Box_Size_Writer, Concrete_Box_Size_Writer, &
    Null_Box_Size_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(writer, periodic_box, box_size_changes)
        class(Abstract_Box_Size_Writer), allocatable, intent(out) :: writer
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: box_size_changes

        if (box_size_changes) then
            allocate(Concrete_Box_Size_Writer :: writer)
        else
            allocate(Null_Box_Size_Writer :: writer)
        end if
        !call writer%construct(periodic_box, )
    end subroutine create

    subroutine destroy(writer)
        class(Abstract_Box_Size_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy

end module procedures_box_size_writer_factory
