module mod_ewald_reci

use, intrinsic :: iso_C_binding, only : C_int, C_double
use data_cell

implicit none

    interface
    
        subroutine C_Epot_reci_nfft_init(C_Ncol) bind(C, name="Epot_reci_nfft_init")
        
            import :: C_int
        
            integer(C_int), value :: C_Ncol
            
        end subroutine C_Epot_reci_nfft_init
        
        subroutine C_Epot_reci_nfft_finalize() bind(C, name="Epot_reci_nfft_finalize")
        
        end subroutine C_Epot_reci_nfft_finalize
        
        subroutine C_Epot_reci_init(C_Lsize, C_alpha) bind(C, name="Epot_reci_init")
        
            import :: Dim, C_double
        
            real(C_double), dimension(Dim) , intent(in) :: C_Lsize
            real(C_double), value :: C_alpha
            
        end subroutine C_Epot_reci_init
        
        function C_Epot_reci_move(C_deltaX, C_Vol) bind(C, name="Epot_reci_move")
        
            import :: Dim, C_double
        
            real(C_double), dimension(Dim), intent(in) :: C_deltaX
            real(C_double), value :: C_Vol
        
        end function C_Epot_reci_move
        
        subroutine C_Epot_reci_updateX(C_lCol, C_xNew) bind(C, name="Epot_reci_updateX")
        
            import :: Dim, C_int, C_double
            
            integer(C_int), value :: C_lCol
            real(C_double), dimension(Dim), intent(in) :: C_xNew
        
        function C_Epot_reci(C_X, C_M, C_Ncol, C_Vol) bind(C, name="Epot_reci")
        
            import :: Dim, C_int, C_double
        
            real(C_double), dimension(Dim, *), intent(in) :: C_X, C_M
            integer(C_int), value :: C_Ncol
            real(C_double), value :: C_Vol
            real(C_double) :: C_Epot_reci
            
        end function C_Epot_reci
    
    end interface

end module mod_ewald_reci
