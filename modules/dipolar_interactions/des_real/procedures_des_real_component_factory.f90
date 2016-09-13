module procedures_des_real_component_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use classes_des_real_component, only: Abstract_DES_Real_Component, Concrete_DES_Real_Component, &
    Null_DES_Real_Component

implicit none

private
public :: create, destroy

contains

    subroutine create(component, periodic_box, box_volume_memento, positions, dipole_moments, &
        interact, pair)
        class(Abstract_DES_Real_Component), allocatable, intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Box_Volume_Memento), intent(in) ::box_volume_memento
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments
        logical, intent(in) :: interact
        class(Abstract_DES_Real_Pair), intent(in) :: pair

        if (interact) then
            allocate(Concrete_DES_Real_Component :: component)
        else
            allocate(Null_DES_Real_Component :: component)
        end if
        call component%construct(periodic_box, box_volume_memento, positions, dipole_moments, pair)
    end subroutine create

    subroutine destroy(component)
        class(Abstract_DES_Real_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy

end module procedures_des_real_component_factory
