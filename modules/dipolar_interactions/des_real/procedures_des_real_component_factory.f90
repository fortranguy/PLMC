module procedures_des_real_component_factory

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use classes_des_real_component, only: Abstract_DES_Real_Component, &
    Concrete_DES_Real_Component, Null_DES_Real_Component
use procedures_property_inquirers, only: components_interact

implicit none

private
public :: create, destroy

contains

    subroutine create(component, periodic_box, positions, dipolar_moments, pair)
        class(Abstract_DES_Real_Component), allocatable, intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        class(Abstract_DES_Real_Pair), intent(in) :: pair

        if (components_interact(pair)) then
            allocate(Concrete_DES_Real_Component :: component)
        else
            allocate(Null_DES_Real_Component :: component)
        end if
        call component%construct(periodic_box, positions, dipolar_moments, pair)
    end subroutine create

    subroutine destroy(component)
        class(Abstract_DES_Real_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy

end module procedures_des_real_component_factory
