module class_ewald_summation_reci

use data_precisions, only: DP
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use module_data, only: test_data_found
use module_types_micro, only: Box_Dimensions, Particle_Index
use module_physics_micro, only: fourier_i
use class_dipolar_hard_spheres

implicit none

private

    type, public :: Ewald_Summation_Reci
    
        integer :: num_wave_vectors
        real(DP), dimension(:, :, :), allocatable :: weight
        complex(DP), dimension(:, :, :), allocatable :: structure
    
    contains
    
        procedure, private :: set_structure => Ewald_Summation_Reci_set_structure
        procedure, private :: set_weight => Ewald_Summation_Reci_set_weight
        !procedure, private :: set_Epot_reci => Ewald_Summation_Reci_set_Epot_reci
        !procedure, private :: Epot_reci_get_structure_modulus => &
                              !Ewald_Summation_Reci_Epot_reci_get_structure_modulus
        !procedure :: reset_Epot_reci_structure => Ewald_Summation_Reci_reset_Epot_reci_structure
        !procedure :: Epot_reci_count_waveVectors => Ewald_Summation_Reci_Epot_reci_count_waveVectors
    
    end type Ewald_Summation_Reci
    
contains

    !> Structure factor init :
    !> \f[
    !>      S(\vec{k}) = \sum_{i} (\vec{k}\cdot\vec{\mu}_i) e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f]
    !> We will also use a restricted definition later :
    !> \f[
    !>      S_\underline{l}(\vec{k}) = \sum_{i \neq l} (\vec{k}\cdot\vec{\mu}_i)
    !>                                 e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f].

    pure subroutine Ewald_Summation_Reci_set_structure(this, Box, this_spheres)

        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3

        real(DP), dimension(Ndim) :: position_div_box, orientation_div_box
        real(DP), dimension(Ndim) :: wave_vector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        integer :: i_particle

        this%structure(:, :, :) = cmplx(0._DP, 0._DP, DP)

        do i_particle = 1, this_spheres%get_num_particles()
        
            position_div_box(:) = 2._DP*PI * this_spheres%get_position(i_particle)/Box%size(:)
            call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
            call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
            call fourier_i(Box%wave(3), position_div_box(3), exp_Ikx_3)
            
            orientation_div_box(:) = this_spheres%get_orientation(i_particle)/Box%size(:)
        
            do kz = -Box%wave(3), Box%wave(3)
                wave_vector(3) = real(kz, DP)

            do ky = -Box%wave(2), Box%wave(2)
                wave_vector(2) = real(ky, DP)
                
            do kx = -Box%wave(1), Box%wave(1)
                wave_vector(1) = real(kx, DP)
                              
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                k_dot_mCol = dot_product(wave_vector, orientation_div_box)
                          
                this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                       cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine Ewald_Summation_Reci_set_structure

    !> \f[
    !>      w(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}{\alpha^2} \sum_{d=1}^3 \frac{k_d^2}{L_d}}}
    !>                                {\sum_{d=1}^3 \frac{k_d^2}{L_d}}
    !> \f]
    
    pure subroutine Ewald_Summation_Reci_set_weight(this, Box, alpha)
        
        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP), intent(in) :: alpha
        
        integer :: kx, ky, kz
        real(DP), dimension(Ndim) :: wave_vector
        real(DP) :: wave_div_box

        do kz = -Box%wave(3), Box%wave(3)
            wave_vector(3) = real(kz, DP)
        
        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = real(ky, DP)
        
        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = real(kx, DP)

            if (kx**2 + ky**2 + kz**2 /= 0) then
                wave_div_box = norm2(wave_vector(:)/Box%size(:))
                this%weight(kx, ky, kz) = exp(-PI**2/alpha**2 * wave_div_box**2) / wave_div_box**2
            else
                this%weight(kx, ky, kz) = 0._DP
            end if

        end do
            
        end do
        
        end do
        
    end subroutine Ewald_Summation_Reci_set_weight

end module class_ewald_summation_reci
