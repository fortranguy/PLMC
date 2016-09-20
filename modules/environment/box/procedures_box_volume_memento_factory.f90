module procedures_box_volume_memento_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento, Retentive_Box_Volume_Memento, &
    Forgetful_Box_Volume_Memento, Null_Box_Volume_Memento
use classes_beta_pressure, only: Abstract_Beta_Pressure
use procedures_environment_inquirers, only: periodicity_is_xyz, box_size_can_change

implicit none

private
public :: create, destroy

contains

    subroutine create(volume_memento, periodic_box, beta_pressure, dipoles_exist)
        class(Abstract_Box_Volume_Memento), allocatable, intent(out) :: volume_memento
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Beta_Pressure), intent(in) :: beta_pressure
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            if (periodicity_is_xyz(periodic_box) .and. box_size_can_change(beta_pressure)) then
                allocate(Retentive_Box_Volume_Memento :: volume_memento)
            else
                allocate(Forgetful_Box_Volume_Memento :: volume_memento)
            end if
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
