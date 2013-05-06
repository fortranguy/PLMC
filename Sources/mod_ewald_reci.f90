module mod_ewald_reci

use, intrinsic :: iso_c_binding, only : c_int, c_double
use data_cell

implicit none

    interface
        
        subroutine c_precalc_exp_ksqr(c_Lsize, c_alpha) bind(c)
        
            import :: Dim, c_double
            real(c_double), intent(in), dimension(Dim) :: c_Lsize
            real(c_double), value :: c_alpha
        
        end subroutine c_precalc_exp_ksqr
        
        subroutine c_nfft_init(Ncol) bind(c)
            import :: c_int
            integer(c_int), value :: Ncol
        end subroutine c_nfft_init
        
        subroutine c_nfft_finalize() bind(c)
            !
        end subroutine c_nfft_finalize
        
        function c_epot_reci(X, D, Ncol, vol) bind(c)
        
            import :: Dim, c_int, c_double
            real(c_double), intent(in), dimension(Dim, *) :: X, D
            integer(c_int), value :: Ncol
            real(c_double), value :: Vol
            real(c_double) :: c_epot_reci
        
        end function c_epot_reci
    
    end interface

end module mod_ewald_reci
