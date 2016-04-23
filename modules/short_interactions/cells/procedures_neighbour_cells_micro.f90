module procedures_neighbour_cells_micro

implicit none

private
public :: pbc_3d_index

contains

    !> Periodic Index, found heuristically.
    pure function pbc_3d_index(ijk_cell, nums_cells)
        integer, intent(in) :: ijk_cell(:), nums_cells(:)
        integer :: pbc_3d_index(3)

        pbc_3d_index = modulo(ijk_cell + nums_cells/2, nums_cells) - nums_cells/2
    end function pbc_3d_index

end module procedures_neighbour_cells_micro
