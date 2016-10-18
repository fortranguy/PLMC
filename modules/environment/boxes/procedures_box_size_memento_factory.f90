module procedures_box_size_memento_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_size_memento, only: Abstract_Box_Size_Memento, Retentive_Box_Size_Memento, &
    Forgetful_Box_Size_Memento, Null_Box_Size_Memento
use procedures_environment_inquirers, only: periodicity_is_xyz

implicit none

private
public :: create, destroy

contains

    subroutine create(memento, periodic_box, box_size_can_change, dipoles_exist)
        class(Abstract_Box_Size_Memento), allocatable, intent(out) :: memento
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: box_size_can_change, dipoles_exist

        if (dipoles_exist) then
            if (periodicity_is_xyz(periodic_box) .and. box_size_can_change) then
                allocate(Retentive_Box_Size_Memento :: memento)
            else
                allocate(Forgetful_Box_Size_Memento :: memento)
            end if
        else
            allocate(Null_Box_Size_Memento :: memento)
        end if
        call memento%construct(periodic_box)
    end subroutine create

    subroutine destroy(memento)
        class(Abstract_Box_Size_Memento), allocatable, intent(inout) :: memento

        if (allocated(memento)) then
            call memento%destroy()
            deallocate(memento)
        end if
    end subroutine destroy

end module procedures_box_size_memento_factory
