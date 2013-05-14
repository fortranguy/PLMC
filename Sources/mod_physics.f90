!> Subroutines for physical purposes

module mod_physics

use data_constants
use data_cell
!$ use omp_lib

implicit none

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    function dist(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2        
        real(DP) :: dist
        
        real(DP), dimension(Dim) :: distVec_12
        
        distVec_12(:) = distVec(xCol1, xCol2)
        
        dist = sqrt(dot_product(distVec_12, distVec_12))
    
    end function dist
    
    function distVec(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP), dimension(Dim) :: distVec
        
        distVec(:) = xCol2(:) - xCol1(:)
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
        
        real(DP) :: gaussSet
        integer :: iset
        
        data iset/0/
        save
        
        if (iset == 0) then
        
            do
            
                call random_number(x)
                x = 2._DP*x - 1._DP
                call random_number(y)
                y = 2._DP*y - 1._DP
            
                normSqr = x*x + y*y
                !write(*, *) "x", x, "y", x, "normSqr", normSqr
                
                if (normSqr <= 1._DP .and. normSqr /= 0._DP) exit
                
            end do
            
            factor = sqrt(-2._DP * log(normSqr) / normSqr)
            
            gaussSet = sigma3d * factor * x
            gauss = sigma3d * factor * y
            
            iset = 1
            
        else
        
            gauss = gaussSet
            iset = 0
            
        end if
        
        !write(*, *) "gauss", gauss
        
    end function gauss
    
    function random_surface()
        
        real(DP), dimension(Dim) :: random_surface
        
        integer :: iDim
        
        do iDim = 1, Dim        
            random_surface(iDim) = gauss()          
        end do
        
        random_surface(:) = random_surface(:) / sqrt(dot_product(random_surface, random_surface))
    
    end function random_surface
    
    subroutine markov_surface(mCol, dm)
    
        real(DP), dimension(Dim), intent(inout) :: mCol
        real(DP), intent(in) :: dm
        
        real(DP), dimension(Dim) :: rotation
        real(DP) :: rotation_dot_mCol
        real(DP) :: amplitude, rand
        integer :: iDim
        
        integer, dimension(2) :: seed

        !call random_seed(get=seed)
        !write(*, *) "markov_seed", seed(:)
        
        do iDim = 1, Dim        
            rotation(iDim) = gauss()
        end do
        
        rotation_dot_mCol = dot_product(rotation, mCol)
        rotation(:) = rotation(:) - rotation_dot_mCol * mCol(:)
        rotation(:) = rotation(:) / sqrt(dot_product(rotation, rotation))
        
        call random_number(rand)
        amplitude = dm * (rand - 0.5_DP)
        
        mCol(:) = mCol(:) + amplitude * rotation(:)
        mCol(:) = mCol(:) / sqrt(dot_product(mCol, mCol))
    
    end subroutine markov_surface

end module mod_physics
