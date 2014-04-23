set_Epot
this%rMin

    subroutine Hard_Spheres_construct_cells(this, Box_size, other, mix_cell_size, mix_range_cut)
        
        class(Hard_Spheres_Potentia), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_range_cut
        
        real(DP), dimension(Ndim) :: same_cell_size
        
        same_cell_size(:) = this%range_cut
        call this%sameCells%construct(Box_size, same_cell_size, this%range_cut)
        call this%sameCells%all_cols_to_cells(this_spheres%get_num_particles(), &
                                              this_spheres%get_all_positions())
        
        call this%mixCells%construct(Box_size, mix_cell_size, mix_range_cut)
        call this%mixCells%all_cols_to_cells(other_spheres%get_num_particles(), &
                                             other_spheres%get_all_positions())

    end subroutine Hard_Spheres_construct_cells

    subroutine Hard_Spheres_Potential_destroy(this)    
        class(Hard_Spheres_Potential), intent(inout) :: this
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()    
    end subroutine Hard_Spheres_Potential_destroy

        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
