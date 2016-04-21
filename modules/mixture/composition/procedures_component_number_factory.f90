module procedures_component_number_factory

use classes_component_number, only: Abstract_Component_Number, Concrete_Component_Number, &
    Null_Component_Number

implicit none

private
public :: create, destroy

contains

    !> Number will be set with coordinates, cf. [[Abstract_Coordinates_Reader]]. Is it too fragile?
    subroutine create(number, exists)
        class(Abstract_Component_Number), allocatable, intent(out) :: number
        logical, intent(in) :: exists

        if (exists) then
            allocate(Concrete_Component_Number :: number)
        else
            allocate(Null_Component_Number :: number)
        end if
    end subroutine create

    subroutine destroy(number)
        class(Abstract_Component_Number), allocatable, intent(inout) :: number

        if (allocated(number)) deallocate(number)
    end subroutine destroy

end module procedures_component_number_factory
