module data_cells

use data_constants, only: num_dimensions

implicit none

private
public :: nums_local_cells

    integer, parameter :: nums_local_cells(num_dimensions) = 3

end module data_cells
