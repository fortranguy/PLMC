module types_cells_wrapper

use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells

implicit none

private

    type, public :: Cells_Wrapper
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
    end type Cells_Wrapper

end module types_cells_wrapper
