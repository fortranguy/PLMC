!> Subroutines for physical purposes

module module_physics

use data_precisions, only : DP
use data_constants, only : PI, sigma3d
use data_box, only : Ndim, Lsize, LsizeMi, Kmax
!$ use omp_lib

implicit none
private
public distVec_PBC, dist_PBC, random_surface, markov_surface, Kmax1_sym, Kmax2_sym, fourier_i, &
       index_from_coord, coord_PBC

contains

    !> Distance between 2 positions with Periodic Boundary Conditions
    
    pure function distVec_PBC(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP), dimension(Ndim) :: distVec_PBC
        
        distVec_PBC(:) = modulo(xCol2(:)-xCol1(:), Lsize(:))
        
        where(distVec_PBC(:) > LsizeMi(:))
            distVec_PBC(:) = distVec_PBC(:) - Lsize(:)
        end where
        
    end function distVec_PBC
    
    pure function dist_PBC(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2        
        real(DP) :: dist_PBC
        
        real(DP), dimension(Ndim) :: distVec_PBC

        distVec_PBC(:) = modulo(xCol2(:)-xCol1(:), Lsize(:))

        where(distVec_PBC(:) > LsizeMi(:))
            distVec_PBC(:) = distVec_PBC(:) - Lsize(:)
        end where
        
        dist_PBC = norm2(distVec_PBC)
    
    end function dist_PBC
    
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
        real(DP) :: amplitude, random
        integer :: iNdim
        
        do iNdim = 1, Ndim        
            rotation(iNdim) = gauss()
        end do
        
        rotation_dot_mCol = dot_product(rotation, mCol)
        rotation(:) = rotation(:) - rotation_dot_mCol * mCol(:)
        rotation(:) = rotation(:) / norm2(rotation)
        
        call random_number(random)
        amplitude = deltaM * (random - 0.5_DP)
        
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
    
    pure subroutine fourier_i(Kmax_i, xColOverL_i, exp_Ikx_i)
    
        integer, intent(in) :: Kmax_i
        real(DP), intent(in) :: xColOverL_i
        complex(DP), dimension(-Kmax_i:Kmax_i), intent(out) :: exp_Ikx_i
        
        integer :: kx_i
        
        exp_Ikx_i(0) = (1._DP, 0._DP)
        exp_Ikx_i(1) = cmplx(cos(xColOverL_i), sin(xColOverL_i), DP)
        exp_Ikx_i(-1) =  conjg(exp_Ikx_i(1))
        
        do kx_i = 2, Kmax_i
        
            exp_Ikx_i(kx_i) = exp_Ikx_i(kx_i-1) * exp_Ikx_i(1)
            exp_Ikx_i(-kx_i) = conjg(exp_Ikx_i(kx_i))
        
        end do
    
    end subroutine fourier_i

    !> 3D index to 1D index
    
    pure function index_from_coord(cell_coord, maxCell_coord)
    
        integer, dimension(:), intent(in) :: cell_coord, maxCell_coord
        integer :: index_from_coord
        
        index_from_coord = cell_coord(1) + maxCell_coord(1)*(cell_coord(2)-1) + &
                           maxCell_coord(1)*maxCell_coord(2)*(cell_coord(3)-1)
    
    end function index_from_coord
    
    !> 3d index Periodic Boundary Conditions
    
    pure function coord_PBC(cell_coord,  maxCell_coord)
    
        integer, dimension(:), intent(in) :: cell_coord, maxCell_coord
        integer, dimension(Ndim) :: coord_PBC
        
        coord_PBC(:) = cell_coord(:)
        
        where (coord_PBC(:) < 1)
            coord_PBC(:) = coord_PBC(:) + maxCell_coord(:)
        end where
        
        where (coord_PBC(:) > maxCell_coord(:))
            coord_PBC(:) = coord_PBC(:) - maxCell_coord(:)
        end where
    
    end function coord_PBC

end module module_physics
