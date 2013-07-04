!> Subroutines for physical purposes

module mod_physics

use data_precisions, only : DP
use data_constants, only : PI, sigma3d
use data_cell, only : Dim, Lsize, LsizeMi, Kmax
!$ use omp_lib

implicit none

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    function dist(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2        
        real(DP) :: dist
        
        real(DP), dimension(Dim) :: distVec_12
        
        distVec_12(:) = distVec(xCol1, xCol2)
        
        dist = norm2(distVec_12)
    
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
        
        real(DP), save :: gaussSave        
        integer, save :: switch = 0
        
        if (switch == 0) then
        
            do
            
                call random_number(x)
                x = 2._DP*x - 1._DP
                call random_number(y)
                y = 2._DP*y - 1._DP
            
                normSqr = x*x + y*y
                
                if (normSqr <= 1._DP .and. normSqr /= 0._DP) exit
                
            end do
            
            factor = sqrt(-2._DP * log(normSqr) / normSqr)
            
            gaussSave = sigma3d * factor * x
            gauss = sigma3d * factor * y
            
            switch = 1
            
        else
        
            gauss = gaussSave
            switch = 0
            
        end if
        
    end function gauss
    
    function random_surface()
        
        real(DP), dimension(Dim) :: random_surface
        
        integer :: iDim
        
        do iDim = 1, Dim        
            random_surface(iDim) = gauss()          
        end do
        
        random_surface(:) = random_surface(:) / norm2(random_surface)
    
    end function random_surface
    
    subroutine markov_surface(mCol, deltaM)
    
        real(DP), dimension(Dim), intent(inout) :: mCol
        real(DP), intent(in) :: deltaM
        
        real(DP), dimension(Dim) :: rotation
        real(DP) :: rotation_dot_mCol
        real(DP) :: amplitude, rand
        integer :: iDim
        
        do iDim = 1, Dim        
            rotation(iDim) = gauss()
        end do
        
        rotation_dot_mCol = dot_product(rotation, mCol)
        rotation(:) = rotation(:) - rotation_dot_mCol * mCol(:)
        rotation(:) = rotation(:) / norm2(rotation)
        
        call random_number(rand)
        amplitude = deltaM * (rand - 0.5_DP)
        
        mCol(:) = mCol(:) + amplitude * rotation(:)
        mCol(:) = mCol(:) / norm2(mCol)
    
    end subroutine markov_surface
    
    !> Fourier coefficients (bases)
    
    subroutine fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
    
        real(DP), dimension(Dim), intent(in) :: xColOverL
        complex(DP), dimension(-Kmax(1):Kmax(1)), intent(out) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)), intent(out) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)), intent(out) :: exp_Ikx_3
        
        real(DP), dimension(Dim) :: arg
        integer :: kx, ky, kz

        arg(:) = 2._DP*PI * 1._DP*xColOverL(:)

        !x
        
        exp_Ikx_1(0) = (1._DP, 0._DP)
        exp_Ikx_1(1) = cmplx(cos(arg(1)), sin(arg(1)), DP)
        exp_Ikx_1(-1) =  conjg(exp_Ikx_1(1))
        
        do kx = 2, Kmax(1)
        
            exp_Ikx_1(kx) = exp_Ikx_1(kx-1) * exp_Ikx_1(1)
            exp_Ikx_1(-kx) = conjg(exp_Ikx_1(kx))
        
        end do
        
        !y
        
        exp_Ikx_2(0) = (1._DP, 0._DP)
        exp_Ikx_2(1) = cmplx(cos(arg(2)), sin(arg(2)), DP)
        exp_Ikx_2(-1) =  conjg(exp_Ikx_2(1))
        
        do ky = 2, Kmax(2)
        
            exp_Ikx_2(ky) = exp_Ikx_2(ky-1) * exp_Ikx_2(1)
            exp_Ikx_2(-ky) = conjg(exp_Ikx_2(ky))
        
        end do
        
        !z
        
        exp_Ikx_3(0) = (1._DP, 0._DP)
        exp_Ikx_3(1) = cmplx(cos(arg(3)), sin(arg(3)), DP)
        exp_Ikx_3(-1) =  conjg(exp_Ikx_3(1))
        
        do kz = 2, Kmax(3)
        
            exp_Ikx_3(kz) = exp_Ikx_3(kz-1) * exp_Ikx_3(1)
            exp_Ikx_3(-kz) = conjg(exp_Ikx_3(kz))
        
        end do        
    
    end subroutine fourier

end module mod_physics
