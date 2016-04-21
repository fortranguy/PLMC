module procedures_component_coordinates_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_number, only: Abstract_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations, Null_Component_Coordinates

implicit none

private
public :: create_positions, create_orientations, destroy

contains

    subroutine create_positions(positions, exists, periodic_box, number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: positions
        logical, intent(in) :: exists
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Number), intent(in) :: number

        if (exists) then
            allocate(Concrete_Component_Positions :: positions)
        else
            allocate(Null_Component_Coordinates :: positions)
        end if
        select type (positions)
            type is (Concrete_Component_Positions)
                call positions%construct(periodic_box, number)
            type is (Null_Component_Coordinates)
                call positions%construct()
            class default
                call error_exit("procedures_component_coordinates_factory: create_positions: "//&
                    "positions: unknown type.")
        end select
    end subroutine create_positions

    subroutine create_orientations(orientations, number, is_dipolar)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: orientations
        class(Abstract_Component_Number), intent(in) :: number
        logical, intent(in) :: is_dipolar

        if (is_dipolar) then
            allocate(Concrete_Component_Orientations :: orientations)
        else
            allocate(Null_Component_Coordinates :: orientations)
        end if
        select type (orientations)
            type is (Concrete_Component_Orientations)
                call orientations%construct(number)
            type is (Null_Component_Coordinates)
                call orientations%destroy()
            class default
                call error_exit("procedures_component_coordinates_factory: create_orientations: "//&
                    "orientations: unknown type.")
        end select
    end subroutine create_orientations

    subroutine destroy(coordinates)
        class(Abstract_Component_Coordinates), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy

end module procedures_component_coordinates_factory
