set_Epot
this%rMin

    subroutine Hard_Spheres_Potential_destroy(this)    
        class(Hard_Spheres_Potential), intent(inout) :: this
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()    
    end subroutine Hard_Spheres_Potential_destroy

        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
