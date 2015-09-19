module procedures_visitable_cells

use data_geometry, only: num_dimensions

implicit none

private
public pbc_3d_index, local_reindex

contains

    pure function pbc_3d_index(i_cell, nums_cells)
        integer, intent(in) :: i_cell(:), nums_cells(:)
        integer :: pbc_3d_index(3)

        pbc_3d_index = modulo(i_cell + nums_cells/2, nums_cells) - nums_cells/2
    end function pbc_3d_index

    pure function local_reindex(i_cell, nums_cells)
        integer, intent(in) :: i_cell(:), nums_cells(:)
        integer :: local_reindex(num_dimensions)

        local_reindex = mod(i_cell, nums_cells) - 1
    end function local_reindex
    ! To find overlap faster

end module procedures_visitable_cells
