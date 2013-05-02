!> Subroutines for physical purposes

module mod_physics

use data_constants
use data_cell
!$ use omp_lib

implicit none

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    function dist_vec(X1, X2)
    
        real(DP), dimension(Dim), intent(in) :: X1, X2

        real(DP), dimension(Dim) :: dist_vec
        
        dist_vec(:) = X2(:) - X1(:)
        dist_vec(:) = modulo(dist_vec(:), Lsize(:))
        
        where( dist_vec(:) > LsizeMi(:) )
            dist_vec(:) = dist_vec(:) - Lsize(:)
        end where
        
    end function dist_vec
    
    function dist(X1, X2)
    
        real(DP), dimension(Dim), intent(in) :: X1, X2
        
        real(DP) :: dist
        real(DP), dimension(Dim) :: DeltaX
        
        DeltaX(:) = dist_vec(X1, X2)
        
        dist = sqrt(dot_product(DeltaX, DeltaX))
    
    end function dist

end module mod_physics
