!> Subroutines for Physics / micro: before

module module_physics_micro

use, intrinsic :: iso_fortran_env, only: error_unit
use data_precisions, only: DP
use data_constants, only: PI, sigma3d
use data_box, only: Ndim
!$ use omp_lib

implicit none
private
public set_discrete_length, sphere_volume, PBC_vector, PBC_distance, random_surface, markov_surface, &
       ewald_real_B, ewald_real_C, &
       NwaveVectors, Box_wave1_sym, Box_wave2_sym, fourier_i, exchange_sign, &
       index_from_coord, coord_PBC, &
       Epot_lennardJones, Epot_yukawa

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
    
    pure function PBC_vector(Box_size, xCol1, xCol2)
    
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP), dimension(Ndim) :: PBC_vector
        
        PBC_vector(:) = modulo(xCol2(:)-xCol1(:), Box_size(:))
        
        where(PBC_vector(:) > Box_size(:)/2._DP)
            PBC_vector(:) = PBC_vector(:) - Box_size(:)
        end where
        
    end function PBC_vector
    
    pure function PBC_distance(Box_size, xCol1, xCol2)
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: xCol1, xCol2
        real(DP) :: PBC_distance
        PBC_distance = norm2(PBC_vector(Box_size, xCol1, xCol2))
    end function PBC_distance
    
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
        
        integer :: jDim
        
        do jDim = 1, Ndim
            random_surface(jDim) = gauss()
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
        integer :: jDim
        
        do jDim = 1, Ndim
            rotation(jDim) = gauss()
        end do
        
        rotation_dot_mCol = dot_product(rotation, mCol)
        rotation(:) = rotation(:) - rotation_dot_mCol * mCol(:)
        rotation(:) = rotation(:) / norm2(rotation)
        
        call random_number(random)
        amplitude = deltaM * (random - 0.5_DP)
        
        mCol(:) = mCol(:) + amplitude * rotation(:)
        mCol(:) = mCol(:) / norm2(mCol)
    
    end subroutine markov_surface
    
    !> \f[ B(r) = \frac{\mathrm{erfc}(\alpha r)}{r^3} +
    !>           2\frac{\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 r^2}}{r^2} \f]
    
    pure function ewald_real_B(alpha, r)    
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r
        real(DP) :: ewald_real_B
    
        ewald_real_B = erfc(alpha*r)/r**3 + 2._DP*alpha/sqrt(PI) * exp(-alpha**2*r**2) / r**2
    
    end function ewald_real_B
    
    !> \f[ C(r) = 3\frac{\mathrm{erfc}(\alpha r)}{r^5} +
    !>            2\frac{\alpha}{\sqrt{\pi}}\left(2\alpha^2 + \frac{3}{r^2}\right)
    !>                                     \frac{e^{-\alpha^2 r^2}}{r^2} \f]
    
    pure function ewald_real_C(alpha, r)    
        real(DP), intent(in) :: alpha
        real(DP), intent(in) :: r
        real(DP) :: ewald_real_C
    
        ewald_real_C = 3._DP*erfc(alpha*r)/r**5 + &
                       2._DP*alpha/sqrt(PI) * (2._DP*alpha**2+3._DP/r**2) * &
                                              exp(-alpha**2*r**2) / r**2
    
    end function ewald_real_C
    
    !> Number of wave vectors
    
    pure function NwaveVectors(Kmax)
    
        integer, dimension(:), intent(in) :: Kmax
        integer :: NwaveVectors
        
        integer :: jDim
        
        NwaveVectors = 1
        do jDim = 1, Ndim
            NwaveVectors = NwaveVectors * (2*Kmax(jDim) + 1)
        end do
    
    end function NwaveVectors

    !> Symmetry: half wave vectors in do loop: Kmax1

    pure function Box_wave1_sym(Kmax, ky, kz)

        integer, dimension(:), intent(in) :: Kmax
        integer, intent(in) :: ky, kz
        integer :: Box_wave1_sym

        if (ky == 0 .and. kz == 0) then
            Box_wave1_sym = 0
        else
            Box_wave1_sym = Kmax(1)
        end if

    end function Box_wave1_sym

    !> Symmetry: half wave vectors in do loop: Kmax2

    pure function Box_wave2_sym(Kmax, kz)
    
        integer, dimension(:), intent(in) :: Kmax
        integer, intent(in) :: kz
        integer :: Box_wave2_sym

        if (kz == 0) then
            Box_wave2_sym = 0
        else
            Box_wave2_sym = Kmax(2)
        end if

    end function Box_wave2_sym
    
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
    
    !> Exchange : + or -
    
    pure function exchange_sign(add)
        
        logical, intent(in) :: add
        real(DP) :: exchange_sign
        
        if (add) then
            exchange_sign = +1._DP
        else
            exchange_sign = -1._DP
        end if
        
    end function exchange_sign

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
