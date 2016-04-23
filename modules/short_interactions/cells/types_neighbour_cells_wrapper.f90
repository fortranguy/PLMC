module types_neighbour_cells_wrapper

use classes_neighbour_cells, only: Abstract_Neighbour_Cells

implicit none

private

    type, public :: Neighbour_Cells_Wrapper
        class(Abstract_Neighbour_Cells), allocatable :: cells
    end type Neighbour_Cells_Wrapper

    type, public :: Neighbour_Cells_Line
        type(Neighbour_Cells_Wrapper), allocatable :: line(:)
    end type Neighbour_Cells_Line

end module types_neighbour_cells_wrapper
