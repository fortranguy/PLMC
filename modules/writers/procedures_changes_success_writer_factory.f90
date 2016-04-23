module procedures_changes_success_writer_factory

use classes_changes_success_writer, only: Concrete_Changes_Selector, &
    Abstract_Changes_Success_Writer, Concrete_Changes_Success_Writer, Null_Changes_Success_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(changes, selector, filename)
        class(Abstract_Changes_Success_Writer), allocatable, intent(out) :: changes
        type(Concrete_Changes_Selector), intent(in) :: selector
        character(len=*), intent(in) :: filename

        if (selector%write_positions) then
            allocate(Concrete_Changes_Success_Writer :: changes)
        else
            allocate(Null_Changes_Success_Writer :: changes)
        end if
        call changes%construct(selector, filename)
    end subroutine create

    subroutine destroy(changes)
        class(Abstract_Changes_Success_Writer), allocatable, intent(inout) :: changes

        if (allocated(changes)) then
            call changes%destroy()
            deallocate(changes)
        end if
    end subroutine destroy

end module procedures_changes_success_writer_factory
