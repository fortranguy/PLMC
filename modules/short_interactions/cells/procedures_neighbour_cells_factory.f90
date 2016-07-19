module procedures_neighbour_cells_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_hard_contact, only: Abstract_Hard_Contact
use classes_pair_potential, only: Abstract_Pair_Potential
use classes_neighbour_cells, only: Abstract_Neighbour_Cells, XYZ_PBC_Neighbour_Cells, &
    XY_PBC_Neighbour_Cells, Null_Neighbour_Cells
use procedures_property_inquirers, only: periodicity_is_xyz, periodicity_is_xy

implicit none

private
public :: create, destroy

contains

    subroutine create(cells, periodic_box, hard_contact, pair_potential, interact)
        class(Abstract_Neighbour_Cells), allocatable, intent(out) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Hard_Contact), intent(in) :: hard_contact
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(in) :: interact

        if (interact) then
            if (periodicity_is_xyz(periodic_box)) then
                allocate(XYZ_PBC_Neighbour_Cells :: cells)
            else if (periodicity_is_xy(periodic_box)) then
                allocate(XY_PBC_Neighbour_Cells :: cells)
            else
                call error_exit("procedures_neighbour_cells_factory: create: "//&
                        "box periodicity is unknown.")
            end if
        else
            allocate(Null_Neighbour_Cells :: cells)
        end if
        call cells%construct(periodic_box, hard_contact, pair_potential)
    end subroutine create

    subroutine destroy(cells)
        class(Abstract_Neighbour_Cells), allocatable, intent(inout) :: cells

        if (allocated(cells)) then
            call cells%destroy()
            deallocate(cells)
        end if
    end subroutine destroy

end module procedures_neighbour_cells_factory
