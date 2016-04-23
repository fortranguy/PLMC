module procedures_neighbour_cells_factory

use procedures_errors, only: error_exit
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use classes_pair_potential, only: Abstract_Pair_Potential
use classes_neighbour_cells, only: Abstract_Neighbour_Cells, XYZ_PBC_Neighbour_Cells, &
    XY_PBC_Neighbour_Cells, Null_Neighbour_Cells

implicit none

private
public :: create, destroy

contains

    subroutine create(cells, periodic_box, pair_potential, interact)
        class(Abstract_Neighbour_Cells), allocatable, intent(out) :: cells
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        logical, intent(in) :: interact

        if (interact) then
            select type (periodic_box)
                type is (XYZ_Periodic_Box)
                    allocate(XYZ_PBC_Neighbour_Cells :: cells)
                type is (XY_Periodic_Box)
                    allocate(XY_PBC_Neighbour_Cells :: cells)
                class default
                    call error_exit("procedures_neighbour_cells_factory: create: "//&
                        "periodic_box type is unknown.")
            end select
        else
            allocate(Null_Neighbour_Cells :: cells)
        end if
        call cells%construct(periodic_box, pair_potential)
    end subroutine create

    subroutine destroy(cells)
        class(Abstract_Neighbour_Cells), allocatable, intent(inout) :: cells

        if (allocated(cells)) then
            call cells%destroy()
            deallocate(cells)
        end if
    end subroutine destroy

end module procedures_neighbour_cells_factory
