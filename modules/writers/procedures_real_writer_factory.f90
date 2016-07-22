module procedures_real_writer_factory

use classes_real_writer, only: Abstract_Real_Writer, Concrete_Real_Writer, Null_Real_Writer

implicit none

private
public :: create, destroy

contains

    subroutine create(writer, needed, filename)
        class(Abstract_Real_Writer), allocatable, intent(out) :: writer
        logical, intent(in) :: needed
        character(len=*), intent(in) :: filename

        if (needed) then
            allocate(Concrete_Real_Writer :: writer)
        else
            allocate(Null_Real_Writer :: writer)
        end if
        call writer%construct(filename)
    end subroutine create

    subroutine destroy(writer)
        class(Abstract_Real_Writer), allocatable, intent(inout) :: writer

        if (allocated(writer)) then
            call writer%destroy()
            deallocate(writer)
        end if
    end subroutine destroy

end module procedures_real_writer_factory
