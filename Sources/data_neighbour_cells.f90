module data_neighbour_cells

use data_box, only: num_dimensions

implicit none

    integer, dimension(num_dimensions), parameter :: num_near_cells_dim = 3
    !< Number of nearest neighbour cells in each direction
    integer, parameter :: num_near_cells = num_near_cells_dim(1) * num_near_cells_dim(2) * &
                                           num_near_cells_dim(3) 
    !< Total number of nearest neighbour cells, including itself

end module data_neighbour_cells
