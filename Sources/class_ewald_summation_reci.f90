module class_ewald_summation_reci

use data_precisions, only: DP
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use module_data, only: test_data_found
use module_types_micro, only: Box_Dimensions, Particle_Index
use module_physics_micro, only: fourier_i, Box_wave1_sym, Box_wave2_sym
use class_dipolar_hard_spheres

implicit none

private

    type, public :: Ewald_Summation_Reci
    
        integer :: num_wave_vectors
        real(DP), dimension(:, :, :), allocatable :: weight
        complex(DP), dimension(:, :, :), allocatable :: structure
    
    contains
    
        procedure :: construct => Ewald_Summation_Reci_construct
        procedure, private :: set_weight => Ewald_Summation_Reci_set_weight
        procedure, private :: set_structure => Ewald_Summation_Reci_set_structure
        procedure :: reset_structure => Ewald_Summation_Reci_reset_structure
        procedure, private :: get_structure_modulus => Ewald_Summation_Reci_get_structure_modulus
        procedure :: count_wave_vectors => Ewald_Summation_Reci_count_wave_vectors
        
        procedure :: move => Ewald_Summation_Reci_move
        procedure :: update_structure_move => Ewald_Summation_Reci_update_structure_move
    
    end type Ewald_Summation_Reci
    
contains

    pure subroutine Ewald_Summation_Reci_construct(this, Box, alpha, this_spheres)
        
        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP), intent(in) :: alpha
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        
        if (allocated(this%weight)) deallocate(this%weight)
        allocate(this%weight(-Box%wave(1):Box%wave(1), &
                             -Box%wave(2):Box%wave(2), &
                             -Box%wave(3):Box%wave(3)))
        
        if (allocated(this%structure)) deallocate(this%structure)
        allocate(this%structure(-Box%wave(1):Box%wave(1), &
                                -Box%wave(2):Box%wave(2), &
                                -Box%wave(3):Box%wave(3)))
        
        call this%set_weight(Box, alpha)
        call this%set_structure(Box, this_spheres)
    
    end subroutine Ewald_Summation_Reci_construct
    
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
        real(DP) :: wave_dot_orientation
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
                wave_dot_orientation = dot_product(wave_vector, orientation_div_box)
                          
                this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                             cmplx(wave_dot_orientation, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine Ewald_Summation_Reci_set_structure
    
    subroutine Ewald_Summation_Reci_reset_structure(this, Box, this_spheres, i_step, modulus_unit)
    
        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: i_step
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reinit
        
        modulus_drifted = this%get_structure_modulus(Box%wave)
        call this%set_structure(Box, this_spheres)
        modulus_reinit = this%get_structure_modulus(Box%wave)
        
        write(modulus_unit, *) i_step, abs(modulus_reinit - modulus_drifted)
    
    end subroutine Ewald_Summation_Reci_reset_structure

    pure function Ewald_Summation_Reci_get_structure_modulus(this, Box_wave) &
                  result(get_structure_modulus)

        class(Ewald_Summation_Reci), intent(in) :: this
        integer, dimension(:), intent(in) :: Box_wave
        real(DP) :: get_structure_modulus

        integer :: kx, ky, kz

        get_structure_modulus = 0._DP

        do kz = 0, Box_wave(3)
            do ky = -Box_wave2_sym(Box_wave, kz), Box_wave(2)
                do kx = -Box_wave1_sym(Box_wave, ky, kz), Box_wave(1)
                    get_structure_modulus = get_structure_modulus + abs(this%structure(kx, ky, kz))
                end do
            end do
        end do

    end function Ewald_Summation_Reci_get_structure_modulus
    
    subroutine Ewald_Summation_Reci_count_wave_vectors(this, Box_wave, wave_unit)

        class(Ewald_Summation_Reci), intent(inout) :: this
        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: wave_unit
        
        integer :: kx, ky, kz

        this%num_wave_vectors = 0

        do kz = 0, Box_wave(3)
            do ky = -Box_wave2_sym(Box_wave, kz), Box_wave(2)
                do kx = -Box_wave1_sym(Box_wave, ky, kz), Box_wave(1)
                    if (kx**2 + ky**2 + kz**2 /= 0) then

                        write(wave_unit, *) kx, ky, kz
                        write(wave_unit, *)
                        write(wave_unit, *)

                        this%num_wave_vectors = this%num_wave_vectors + 1

                    end if
                end do
            end do
        end do

    end subroutine Ewald_Summation_Reci_count_wave_vectors
    
    !> Move

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta S^2
    !> w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta S^2 = 2\Re[
    !>                  (\vec{\mu}_l\cdot\vec{k})
    !>                  (e^{-i\vec{k}\cdot\vec{x}^\prime_l} - e^{-i\vec{k}\cdot\vec{x}_l})
    !>                  S_\underline{l}(\vec{k})
    !>               ]
    !> \f]

    pure function Ewald_Summation_Reci_move(this, Box, old, new) &
                  result(move)

        class(Ewald_Summation_Reci), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: move
        
        real(DP), dimension(Ndim) :: new_position_div_box, old_position_div_box
        real(DP), dimension(Ndim) :: orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP) :: real_part
        
        real(DP), dimension(Ndim) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)
        
        old_position_div_box(:) = 2._DP*PI * old%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        orientation_div_box(:) = new%orientation(:)/Box%size(:)

        move = 0._DP

        do kz = 0, Box%wave(3) ! symmetry: half wave vectors -> double Energy
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)
                    
                    wave_dot_orientation = dot_product(wave_vector, orientation_div_box)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)

                    real_part = wave_dot_orientation * real((conjg(exp_IkxNew) - conjg(exp_IkxOld)) * &
                                (this%structure(kx, ky, kz) - cmplx(wave_dot_orientation, 0._DP, DP) * &
                                exp_IkxOld), DP)

                    move = move + 2._DP * this%weight(kx, ky, kz) * real_part

                end do
            
            end do
        
        end do

        move = 4._DP*PI/product(Box%size) * move

    end function Ewald_Summation_Reci_move

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    pure subroutine Ewald_Summation_Reci_update_structure_move(this, Box, old, new)

        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        
        real(DP), dimension(Ndim) :: new_position_div_box, old_position_div_box
        real(DP), dimension(Ndim) :: orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(Ndim) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)
        
        old_position_div_box(:) = 2._DP*PI * old%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        orientation_div_box(:) = new%orientation(:)/Box%size(:)

        do kz = 0, Box%wave(3)
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)

                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)
                    
                    wave_dot_orientation = dot_product(wave_vector, orientation_div_box)
                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                                                          
                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(wave_dot_orientation, 0._DP, DP) * &
                                                 (exp_IkxNew - exp_IkxOld)

                end do
                
            end do
            
        end do

    end subroutine Ewald_Summation_Reci_update_structure_move

end module class_ewald_summation_reci
