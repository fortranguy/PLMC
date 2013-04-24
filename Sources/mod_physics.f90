!> Subroutines for physical purposes

module mod_physics

use data_constants
use data_cell
!$ use omp_lib

implicit none

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    function dist(X1, X2)
    
        real(DP), dimension(Dim), intent(in) :: X1, X2
        
        real(DP), dimension(Dim) :: DeltaX
        real(DP) :: dist
        
        DeltaX(:) = X2(:) - X1(:)
        DeltaX(:) = modulo(DeltaX(:), Lsize(:))
        
        where( DeltaX(:) > LsizeMi(:) )
            DeltaX(:) = DeltaX(:) - Lsize(:)
        end where
        
        dist = sqrt(dot_product(DeltaX, DeltaX))
    
    end function dist

end module mod_physics
