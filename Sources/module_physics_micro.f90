!> Subroutines for Physics / micro: before

module module_physics_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64, error_unit
use data_constants, only: PI, sigma3d
use data_box, only: num_dimensions
use module_geometry, only: geometry

implicit none
private
public set_discrete_length, sphere_volume, PBC_vector, PBC_distance, random_surface, markov_surface, &
       dipolar_pair_energy, ewald_real_B, ewald_real_C, &
       num_wave_vectors, Box_wave1_sym, Box_wave2_sym, fourier_i, set_exp_kz, exchange_sign, &
       index_from_coord, coord_PBC, &
       potential_energy_lennard_jones, potential_energy_yukawa

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
    
    pure function PBC_vector(Box_size, position1, position2)
    
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: position1, position2
        real(DP), dimension(num_dimensions) :: PBC_vector
        
        if (geometry%bulk) then
        
            PBC_vector(:) = modulo(position2(:)-position1(:), Box_size(:))
            
            where(PBC_vector(:) > Box_size(:)/2._DP)
                PBC_vector(:) = PBC_vector(:) - Box_size(:)
            end where
            
        else if (geometry%slab) then
        
            PBC_vector(1:num_dimensions-1) = &
                modulo(position2(1:num_dimensions-1) - position1(1:num_dimensions-1), &
                       Box_size(1:num_dimensions-1))

            where(PBC_vector(1:num_dimensions-1) > Box_size(1:num_dimensions-1)/2._DP)
                PBC_vector(1:num_dimensions-1) = &
                    PBC_vector(1:num_dimensions-1) - Box_size(1:num_dimensions-1)
            end where

            PBC_vector(num_dimensions) = position2(num_dimensions) - position1(num_dimensions)
        
        end if
        
    end function PBC_vector
    
    pure function PBC_distance(Box_size, position1, position2)
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: position1, position2
        real(DP) :: PBC_distance
        PBC_distance = norm2(PBC_vector(Box_size, position1, position2))
    end function PBC_distance
    
    !> Rotation
    !> From SMAC, Algorithm 1.19, p.39
    
    function gauss()
        
        real(DP) :: gauss
        
        real(DP) :: x, y
        real(DP) :: normSqr, factor
        
        real(DP), save :: gauss_save
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
            
            gauss_save = sigma3d * factor * x
            gauss = sigma3d * factor * y
            
            switch = 1
            
        else
        
            gauss = gauss_save
            switch = 0
            
        end if
        
    end function gauss
    
    !> From SMAC, Algorithm 1.23, p. 43
    
    function random_surface()
        
        real(DP), dimension(num_dimensions) :: random_surface
        
        integer :: i_dim
        
        do i_dim = 1, num_dimensions
            random_surface(i_dim) = gauss()
        end do
        
        random_surface(:) = random_surface(:) / norm2(random_surface)
    
    end function random_surface
    
    !> From SMAC, Algorithm 1.24, p.44
    
    subroutine markov_surface(orientation, orientation_delta)
    
        real(DP), dimension(:), intent(inout) :: orientation
        real(DP), intent(in) :: orientation_delta
        
        real(DP), dimension(num_dimensions) :: rotation
        real(DP) :: rotation_dot_orientation
        real(DP) :: amplitude, random
        integer :: i_dim
        
        do i_dim = 1, num_dimensions
            rotation(i_dim) = gauss()
        end do
        
        rotation_dot_orientation = dot_product(rotation, orientation)
        rotation(:) = rotation(:) - rotation_dot_orientation * orientation(:)
        rotation(:) = rotation(:) / norm2(rotation)
        
        call random_number(random)
        amplitude = orientation_delta * (random - 0.5_DP)
        
        orientation(:) = orientation(:) + amplitude * rotation(:)
        orientation(:) = orientation(:) / norm2(orientation)
    
    end subroutine markov_surface

    !> \f[
    !>      \frac{(\vec{\mu}_i\cdot\vec{\mu_j})}{r^3} - 
    !>     3\frac{(\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij})}{r^5}
    !> \f]
    
    pure function dipolar_pair_energy(orientation_i, orientation_j, vector_ij)
    
        real(DP), dimension(:), intent(in) :: orientation_i, orientation_j
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP) :: dipolar_pair_energy
        
        real(DP) :: distance_ij
        
        distance_ij = norm2(vector_ij)        
        dipolar_pair_energy = dot_product(orientation_i, orientation_j) / distance_ij**3 - &
                              3._DP * dot_product(orientation_i, vector_ij) * &
                                      dot_product(orientation_j, vector_ij) / distance_ij**5
        
    end function dipolar_pair_energy
    
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
    
    pure function num_wave_vectors(Box_wave)
    
        integer, dimension(:), intent(in) :: Box_wave
        integer :: num_wave_vectors
        
        integer :: i_dim
        
        num_wave_vectors = 1
        do i_dim = 1, num_dimensions
            num_wave_vectors = num_wave_vectors * (2*Box_wave(i_dim) + 1)
        end do
    
    end function num_wave_vectors

    !> Symmetry: half wave vectors in do loop: Box_wave1

    pure function Box_wave1_sym(Box_wave, ky, kz)

        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: ky, kz
        integer :: Box_wave1_sym

        if (ky == 0 .and. kz == 0) then
            Box_wave1_sym = 0
        else
            Box_wave1_sym = Box_wave(1)
        end if

    end function Box_wave1_sym

    !> Symmetry: half wave vectors in do loop: Box_wave2

    pure function Box_wave2_sym(Box_wave, kz)
    
        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: kz
        integer :: Box_wave2_sym

        if (kz == 0) then
            Box_wave2_sym = 0
        else
            Box_wave2_sym = Box_wave(2)
        end if

    end function Box_wave2_sym
    
    !> Fourier coefficients (bases)
    
    pure subroutine fourier_i(Box_wave_i, position_div_box_i, exp_Ikx_i)
    
        integer, intent(in) :: Box_wave_i
        real(DP), intent(in) :: position_div_box_i
        complex(DP), dimension(-Box_wave_i:Box_wave_i), intent(out) :: exp_Ikx_i
        
        integer :: kx_i
        
        exp_Ikx_i(0) = (1._DP, 0._DP)
        exp_Ikx_i(1) = cmplx(cos(position_div_box_i), sin(position_div_box_i), DP)
        exp_Ikx_i(-1) = conjg(exp_Ikx_i(1))
        
        do kx_i = 2, Box_wave_i
            exp_Ikx_i(kx_i) = exp_Ikx_i(kx_i-1) * exp_Ikx_i(1)
            exp_Ikx_i(-kx_i) = conjg(exp_Ikx_i(kx_i))
        end do
    
    end subroutine fourier_i
    
    pure subroutine set_exp_kz(Box_wave, wave_norm, zCol, exp_kzCol_tab)
    
        integer, dimension(:), intent(in) :: Box_wave
        real(DP), dimension(-Box_wave(1):Box_wave(1), -Box_wave(2):Box_wave(2)), intent(in) :: wave_norm
        real(DP), intent(in) :: zCol
        real(DP), dimension(0:Box_wave(1), 0:Box_wave(2)), intent(out) :: exp_kzCol_tab

        integer :: kx, ky

        do ky = 0, Box_wave(2)
        
            do kx = ky, Box_wave(1)
                exp_kzCol_tab(kx, ky) = exp(wave_norm(kx, ky) * zCol)
            end do
            
            do kx = 0, ky-1
                exp_kzCol_tab(kx, ky) = exp_kzCol_tab(ky, kx)
            end do
            
        end do

    end subroutine set_exp_kz
    
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
    
    pure function index_from_coord(cell_coord, max_cell_coord)
    
        integer, dimension(:), intent(in) :: cell_coord, max_cell_coord
        integer :: index_from_coord
        
        index_from_coord = cell_coord(1) + max_cell_coord(1)*(cell_coord(2)-1) + &
                           max_cell_coord(1)*max_cell_coord(2)*(cell_coord(3)-1)
    
    end function index_from_coord
    
    !> 3d index Periodic Boundary Conditions
    
    pure function coord_PBC(cell_coord,  max_cell_coord)
    
        integer, dimension(:), intent(in) :: cell_coord, max_cell_coord
        integer, dimension(num_dimensions) :: coord_PBC
        
        coord_PBC(:) = cell_coord(:)
        
        where (coord_PBC(:) < 1)
            coord_PBC(:) = coord_PBC(:) + max_cell_coord(:)
        end where
        
        where (coord_PBC(:) > max_cell_coord(:))
            coord_PBC(:) = coord_PBC(:) - max_cell_coord(:)
        end where
    
    end function coord_PBC
    
    !> Potential energy
    
    !> Lennard-Jones
    !> \f[
    !>      4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
    !> \f]
    pure function potential_energy_lennard_jones(epsilon, sigma, r)
    
        real(DP), intent(in) :: epsilon
        real(DP), intent(in) :: sigma, r
        real(DP) :: potential_energy_lennard_jones
        
        potential_energy_lennard_jones = 4._DP * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    
    end function potential_energy_lennard_jones
    
    !> Yukawa potential_energy
    !> \f[ \epsilon \frac{e^{-\alpha (r-r_0)}}{r} \f]
    
    pure function potential_energy_yukawa(epsilon, alpha, r_0, r)
    
        real(DP), intent(in) :: epsilon, alpha
        real(DP), intent(in) :: r_0, r
        real(DP) :: potential_energy_yukawa
        
        potential_energy_yukawa = epsilon * exp(-alpha*(r-r_0)) / r
    
    end function potential_energy_yukawa

end module module_physics_micro
