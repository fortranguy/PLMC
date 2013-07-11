!> Subroutines for physical purposes

module mod_physics

use data_precisions, only : DP
use data_constants, only : PI, sigma3d
use data_cell, only : Ndim, Lsize, LsizeMi, Kmax
!$ use omp_lib

implicit none
private
public dist, distVec, random_surface, markov_surface, Kmax1_sym, Kmax2_sym, fourier, index_from_coord

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    pure function dist(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2        
        real(DP) :: dist
        
        real(DP), dimension(Ndim) :: distVec_12
        
        distVec_12(:) = distVec(xCol1, xCol2)
        
        dist = norm2(distVec_12)
    
    end function dist
    
    pure function distVec(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP), dimension(Ndim) :: distVec
        
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
        
        real(DP), dimension(Ndim) :: random_surface
        
        integer :: iNdim
        
        do iNdim = 1, Ndim        
            random_surface(iNdim) = gauss()          
        end do
        
        random_surface(:) = random_surface(:) / norm2(random_surface)
    
    end function random_surface
    
    subroutine markov_surface(mCol, deltaM)
    
        real(DP), dimension(:), intent(inout) :: mCol
        real(DP), intent(in) :: deltaM
        
        real(DP), dimension(Ndim) :: rotation
        real(DP) :: rotation_dot_mCol
        real(DP) :: amplitude, rand
        integer :: iNdim
        
        do iNdim = 1, Ndim        
            rotation(iNdim) = gauss()
        end do
        
        rotation_dot_mCol = dot_product(rotation, mCol)
        rotation(:) = rotation(:) - rotation_dot_mCol * mCol(:)
        rotation(:) = rotation(:) / norm2(rotation)
        
        call random_number(rand)
        amplitude = deltaM * (rand - 0.5_DP)
        
        mCol(:) = mCol(:) + amplitude * rotation(:)
        mCol(:) = mCol(:) / norm2(mCol)
    
    end subroutine markov_surface

    !> Symmetry : half wave vectors in do loop : Kmax1

    pure function Kmax1_sym(ky, kz)

        integer, intent(in) :: ky, kz
        integer :: Kmax1_sym

        if (ky == 0 .and. kz == 0) then
            Kmax1_sym = 0
        else
            Kmax1_sym = Kmax(1)
        end if

    end function Kmax1_sym

    !> Symmetry : half wave vectors in do loop : Kmax2

    pure function Kmax2_sym(kz)

        integer, intent(in) :: kz
        integer :: Kmax2_sym

        if (kz == 0) then
            Kmax2_sym = 0
        else
            Kmax2_sym = Kmax(2)
        end if

    end function Kmax2_sym
    
    !> Fourier coefficients (bases)
    
    pure subroutine fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
    
        real(DP), dimension(:), intent(in) :: xColOverL
        complex(DP), dimension(-Kmax(1):Kmax(1)), intent(out) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)), intent(out) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)), intent(out) :: exp_Ikx_3
        
        real(DP), dimension(Ndim) :: arg
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
    
    !> 3d vector to 1D index 
    
    pure function index_from_coord(cell_coord, maxCell_coord)
    
        integer, dimension(:), intent(in) :: cell_coord, maxCell_coord
        integer :: index_from_coord
        
        index_from_coord = cell_coord(1) + maxCell_coord(1)*(cell_coord(2)-1) + &
                           maxCell_coord(1)*maxCell_coord(2)*(cell_coord(3)-1)
    
    end function index_from_coord

end module mod_physics
