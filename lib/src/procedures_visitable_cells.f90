module procedures_visitable_cells

implicit none

private
public pbc_3d_index

contains

    pure function pbc_3d_index(i_cell, nums_cells)
        integer, intent(in) :: i_cell(:), nums_cells(:)
        integer :: pbc_3d_index(3)

        pbc_3d_index = modulo(i_cell + nums_cells/2, nums_cells) - nums_cells/2
    end function pbc_3d_index

end module procedures_visitable_cells
