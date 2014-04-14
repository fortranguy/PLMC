

    
    procedure :: adapt_move_delta => Hard_Spheres_adapt_move_delta      
    procedure :: set_move_delta => Hard_Spheres_set_move_delta  

    subroutine Hard_Spheres_adapt_move_delta(this, Box_size, reject)
        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject        
    
        call this%move%adapt_delta(Box_size, reject)
    end subroutine Hard_Spheres_adapt_move_delta
    
    subroutine Hard_Spheres_set_move_delta(this, Box_size, reject, report_unit)
        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        call this%move%set_delta(this%name, Box_size, reject, report_unit)
    end subroutine Hard_Spheres_set_move_delta
