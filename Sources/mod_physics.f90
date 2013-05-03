!> Subroutines for physical purposes

module mod_physics

use data_constants
use data_cell
!$ use omp_lib

implicit none

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    function dist(X1, X2)
    
        real(DP), dimension(:), intent(in) :: X1, X2        
        real(DP) :: dist
        
        real(DP), dimension(Dim) :: DeltaX
        
        DeltaX(:) = distVec(X1, X2)
        
        dist = sqrt(dot_product(DeltaX, DeltaX))
    
    end function dist
    
    function distVec(X1, X2)
    
        real(DP), dimension(:), intent(in) :: X1, X2
        real(DP), dimension(Dim) :: distVec
        
        distVec(:) = X2(:) - X1(:)
        distVec(:) = modulo(distVec(:), Lsize(:))
        
        where(distVec(:) > LsizeMi(:))
            distVec(:) = distVec(:) - Lsize(:)
        end where
        
    end function distVec
    
    !> Rotation
    
    function gauss()
        
        real(DP) :: gauss
        
        real(DP) :: x, y
        real(DP) :: normSqr, factor
        
        do
        
            call random_number(x)
            x = 2._DP*x - 1._DP
            call random_number(y)
            y = 2._DP*y - 1._DP
        
            normSqr = x*x + y*y
            if (normSqr <= 1._DP .and. normSqr /= 0._DP) exit
            
        end do
        
        factor = sqrt(-2._DP * log(normSqr) / normSqr)
        gauss = sigma3d * factor * x
        
    end function gauss
    
    function random_surface()
        
        real(DP), dimension(Dim) :: random_surface
        
        integer :: iDim
        
        do iDim = 1, Dim        
            random_surface(iDim) = gauss()          
        end do
        
        random_surface(:) = random_surface(:) / sqrt(dot_product(random_surface, random_surface))
    
    end function random_surface

end module mod_physics
