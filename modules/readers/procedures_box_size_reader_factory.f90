module procedures_box_size_reader_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_reader, only: Abstract_Box_Size_Reader, Concrete_Box_Size_Reader

implicit none

private
public :: create, destroy

contains

    subroutine create(box_size, periodic_box)
        class(Abstract_Box_Size_Reader), allocatable, intent(out) :: box_size
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        allocate(Concrete_Box_Size_Reader :: box_size)
        call box_size%construct(periodic_box)
    end subroutine create

    subroutine destroy(box_size)
        class(Abstract_Box_Size_Reader), allocatable, intent(inout) :: box_size

        if (allocated(box_size)) then
            call box_size%destroy()
            deallocate(box_size)
        end if
    end subroutine destroy

end module procedures_box_size_reader_factory
