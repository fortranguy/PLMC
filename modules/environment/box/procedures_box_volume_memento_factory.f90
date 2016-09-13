module procedures_box_volume_memento_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento, Concrete_Box_Volume_Memento, &
    Null_Box_Volume_Memento

implicit none

private
public :: create, destroy

contains

    subroutine create(volume_memento, periodic_box, needed)
        class(Abstract_Box_Volume_Memento), allocatable, intent(out) :: volume_memento
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Box_Volume_Memento :: volume_memento)
        else
            allocate(Null_Box_Volume_Memento :: volume_memento)
        end if

        call volume_memento%construct(periodic_box)
    end subroutine create

    subroutine destroy(volume_memento)
        class(Abstract_Box_Volume_Memento), allocatable, intent(inout) :: volume_memento

        if (allocated(volume_memento)) then
            call volume_memento%destroy()
            deallocate(volume_memento)
        end if
    end subroutine destroy

end module procedures_box_volume_memento_factory
