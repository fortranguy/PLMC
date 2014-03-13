!> Subroutines for Physics / micro: before

module module_physics_micro

use, intrinsic :: iso_fortran_env, only: error_unit
use data_precisions, only: DP
use data_constants, only: PI, sigma3d
use data_box, only: Ndim, Lsize, Kmax
!$ use omp_lib

implicit none
private
public set_discrete_length, sphere_volume, distVec_PBC, dist_PBC, random_surface, markov_surface, &
       Kmax1_sym, Kmax2_sym, fourier_i, index_from_coord, coord_PBC, Epot_lennardJones, Epot_yukawa

contains

    !> Set discretization length
    
    subroutine set_discrete_length(lMin, dl)
        
        real(DP), intent(in) ::lMin
        real(DP), intent(inout) :: dl
        
        if (dl > lMin) then
            write(error_unit, *) "    dl", dl, "> lMin", lMin
            dl = lMin
            write(error_unit, *) "    dl <- lMin"
        end if
    
    end subroutine set_discrete_length
    
    !> Calculate the volume of the sphere
    
    pure function sphere_volume(radius)
        real(DP), intent(in) :: radius
        real(DP) :: sphere_volume
        sphere_volume = 4._DP/3._DP * PI * radius**3
    end function sphere_volume

    !> Distance between 2 positions with Periodic Boundary Conditions
    !> from SMAC, algorithm 2.5 & 2.6, p.91
    
    pure function distVec_PBC(xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP), dimension(Ndim) :: distVec_PBC
        
        distVec_PBC(:) = modulo(xCol2(:)-xCol1(:), Lsize(:))
        
        where(distVec_PBC(:) > Lsize(:)/2._DP)
            distVec_PBC(:) = distVec_PBC(:) - Lsize(:)
        end where
        
    end function distVec_PBC
    
    pure function dist_PBC(xCol1, xCol2)
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP) :: dist_PBC
        dist_PBC = norm2(distVec_PBC(xCol1, xCol2))
    end function dist_PBC
    
    !> Rotation
    !> From SMAC, Algorithm 1.19, p.39
    
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
                
                if (normSqr <= 1._DP .and. normSqr > 0._DP) exit
                
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
    
    !> From SMAC, Algorithm 1.23, p. 43
    
    function random_surface()
        
        real(DP), dimension(Ndim) :: random_surface
        
        integer :: iNdim
        
        do iNdim = 1, Ndim
            random_surface(iNdim) = gauss()
        end do
        
        random_surface(:) = random_surface(:) / norm2(random_surface)
    
    end function random_surface
    
    !> From SMAC, Algorithm 1.24, p.44
    
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

    !> Symmetry: half wave vectors in do loop: Kmax1

    pure function Kmax1_sym(ky, kz)

        integer, intent(in) :: ky, kz
        integer :: Kmax1_sym

        if (ky == 0 .and. kz == 0) then
            Kmax1_sym = 0
        else
            Kmax1_sym = Kmax(1)
        end if

    end function Kmax1_sym

    !> Symmetry: half wave vectors in do loop: Kmax2

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
        exp_Ikx_i(-1) = conjg(exp_Ikx_i(1))
        
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
    
    !> Potential energy
    
    !> Lennard-Jones
    !> \f[
    !>      4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
    !> \f]
    pure function Epot_lennardJones(epsilon, sigma, r)
    
        real(DP), intent(in) :: epsilon
        real(DP), intent(in) :: sigma, r
        real(DP) :: Epot_lennardJones
        
        Epot_lennardJones = 4._DP * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    
    end function Epot_lennardJones
    
    !> Yukawa potential
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_0)}}{r} \f]
    
    pure function Epot_yukawa(epsilon, alpha, r_0, r)
    
        real(DP), intent(in) :: epsilon, alpha
        real(DP), intent(in) :: r_0, r
        real(DP) :: Epot_yukawa
        
        Epot_yukawa = epsilon * exp(-alpha*(r-r_0)) / r
    
    end function Epot_yukawa

end module module_physics_micro