module class_ewald_summation_reci

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use module_types_micro, only: Box_Parameters, Particle_Index
use module_physics_micro, only: fourier_i, Box_wave1_sym, Box_wave2_sym, exchange_sign
use class_hard_spheres, only: Dipolar_Hard_Spheres

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
        procedure :: destroy => Ewald_Summation_Reci_destroy
        procedure :: reset_structure => Ewald_Summation_Reci_reset_structure
        procedure, private :: get_structure_modulus => Ewald_Summation_Reci_get_structure_modulus
        procedure :: count_wave_vectors => Ewald_Summation_Reci_count_wave_vectors
        procedure :: total_energy => Ewald_Summation_Reci_total_energy
        
        procedure :: move_energy => Ewald_Summation_Reci_move_energy
        procedure :: update_structure_move => Ewald_Summation_Reci_update_structure_move
        procedure :: rotation_energy => Ewald_Summation_Reci_rotation_energy
        procedure :: update_structure_rotation => Ewald_Summation_Reci_update_structure_rotation
        procedure :: exchange_energy => Ewald_Summation_Reci_exchange_energy
        procedure :: update_structure_exchange => Ewald_Summation_Reci_update_structure_exchange
    
    end type Ewald_Summation_Reci
    
contains

    pure subroutine Ewald_Summation_Reci_construct(this, Box, alpha, this_spheres)
        
        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
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
        type(Box_Parameters), intent(in) :: Box
        real(DP), intent(in) :: alpha
        
        integer :: kx, ky, kz
        real(DP), dimension(num_dimensions) :: wave_vector
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
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres

        complex(DP) :: exp_Ikx
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3

        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: orientation_div_box
        real(DP), dimension(num_dimensions) :: wave_vector
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
                              
                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                wave_dot_orientation = dot_product(wave_vector, orientation_div_box)
                          
                this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                             cmplx(wave_dot_orientation, 0._DP, DP) * exp_Ikx
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine Ewald_Summation_Reci_set_structure
    
    subroutine Ewald_Summation_Reci_destroy(this)
    
        class(Ewald_Summation_Reci), intent(inout) :: this
        
        if (allocated(this%structure)) deallocate(this%structure)
        if (allocated(this%weight)) deallocate(this%weight)
    
    end subroutine Ewald_Summation_Reci_destroy
    
    subroutine Ewald_Summation_Reci_reset_structure(this, Box, this_spheres, i_step, modulus_unit)
    
        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: i_step
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reset
        
        modulus_drifted = this%get_structure_modulus(Box%wave)
        call this%set_structure(Box, this_spheres)
        modulus_reset = this%get_structure_modulus(Box%wave)
        
        write(modulus_unit, *) i_step, abs(modulus_reset - modulus_drifted)
    
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
    
    pure function Ewald_Summation_Reci_total_energy(this, Box) result(total_energy)
        
        class(Ewald_Summation_Reci), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP) :: total_energy

        integer :: kx, ky, kz

        total_energy = 0._DP

        do kz = -Box%wave(3), Box%wave(3)
            do ky = -Box%wave(2), Box%wave(2)
                do kx = -Box%wave(1), Box%wave(1)
                    total_energy = total_energy + this%weight(kx, ky, kz) * &
                                                  real(this%structure(kx, ky, kz) * &
                                                  conjg(this%structure(kx, ky, kz)), DP)
                end do
            end do
        end do
        
        total_energy = 2._DP*PI / product(Box%size) * total_energy
        
    end function Ewald_Summation_Reci_total_energy
    
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

    pure function Ewald_Summation_Reci_move_energy(this, Box, old, new) result(move_energy)

        class(Ewald_Summation_Reci), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: move_energy
        
        real(DP), dimension(num_dimensions) :: new_position_div_box, old_position_div_box
        real(DP), dimension(num_dimensions) :: orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP) :: real_part
        
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)
        
        old_position_div_box(:) = 2._DP*PI * old%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        orientation_div_box(:) = new%orientation(:) / Box%size(:)

        move_energy = 0._DP

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

                    move_energy = move_energy + 2._DP * this%weight(kx, ky, kz) * real_part

                end do
            
            end do
        
        end do

        move_energy = 4._DP*PI / product(Box%size) * move_energy

    end function Ewald_Summation_Reci_move_energy

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    pure subroutine Ewald_Summation_Reci_update_structure_move(this, Box, old, new)

        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        
        real(DP), dimension(num_dimensions) :: new_position_div_box, old_position_div_box
        real(DP), dimension(num_dimensions) :: orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)
        
        old_position_div_box(:) = 2._DP*PI * old%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        orientation_div_box(:) = new%orientation(:) / Box%size(:)

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
    
    !> Rotation

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta S^2
    !>                                       w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta S^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2\Re\{
    !>                  [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>                  e^{-i \vec{k} \cdot \vec{x}_l} S_\underline{l}(\vec{k})
    !>               \}
    !> \f]
    
    pure function Ewald_Summation_Reci_rotation_energy(this, Box, old, new) result(rotation_energy)

        class(Ewald_Summation_Reci), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation_energy

        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: new_orientation_div_box, old_orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP) :: real_part
        
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_new_orientation, wave_dot_old_orientation
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call fourier_i(Box%wave(3), position_div_box(3), exp_Ikx_3)

        new_orientation_div_box(:) = new%orientation(:) / Box%size(:)
        old_orientation_div_box(:) = old%orientation(:) / Box%size(:)

        rotation_energy = 0._DP

        do kz = 0, Box%wave(3)
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)

                    wave_dot_new_orientation = dot_product(wave_vector, new_orientation_div_box)
                    wave_dot_old_orientation = dot_product(wave_vector, old_orientation_div_box)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    real_part = wave_dot_new_orientation**2 - wave_dot_old_orientation**2
                    real_part = real_part + &
                                2._DP * (wave_dot_new_orientation - wave_dot_old_orientation) * &
                                real(conjg(exp_Ikx) * &
                                (this%structure(kx, ky, kz) - &
                                wave_dot_old_orientation * exp_Ikx), DP)

                    rotation_energy = rotation_energy + this%weight(kx, ky, kz) * real_part

                end do

            end do

        end do

        rotation_energy = 4._DP*PI / product(Box%size) * rotation_energy

    end function Ewald_Summation_Reci_rotation_energy

    !> Update moment -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = [\vec{k}\cdot(\vec{\mu}_l^\prime - \vec{\mu}_l)]
    !>                          e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    pure subroutine Ewald_Summation_Reci_update_structure_rotation(this, Box, old, new)

        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new

        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: new_orientation_div_box, old_orientation_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: k_dot_deltaMcol
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call fourier_i(Box%wave(3), position_div_box(3), exp_Ikx_3)

        new_orientation_div_box(:) = new%orientation(:)/Box%size(:)
        old_orientation_div_box(:) = old%orientation(:)/Box%size(:)

        do kz = 0, Box%wave(3)
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)

                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)
                    
                    k_dot_deltaMcol = dot_product(wave_vector, &
                                      new_orientation_div_box - old_orientation_div_box)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_Ikx

                end do
                
            end do
            
        end do

    end subroutine Ewald_Summation_Reci_update_structure_rotation
    
    !> Energy of 1 dipole with others
    
    !> Addition :
    
    !> Difference of Energy
    !> \f[ \Delta U^{N+1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}}
    !>                          (\vec{k} \cdot +\vec{\mu}_{N+1}) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot +\vec{\mu}_{N+1}) +
    !>                              2\Re[S(\vec{k}) e^{-i \vec{k} \cdot \vec{x}_{N+1}}]
    !>                          \}
    !> \f]
    
    !> Summary: only the sign of \f$\vec{\mu}\f$ changes.

    pure function Ewald_Summation_Reci_exchange_energy(this, Box, particle) result(exchange_energy)

        class(Ewald_Summation_Reci), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange_energy
        
        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: orientation_div_box
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx
        
        real(DP) :: real_part
        
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz
        
        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call fourier_i(Box%wave(3), position_div_box(3), exp_Ikx_3)
        
        orientation_div_box(:) = exchange_sign(particle%add) * particle%orientation(:) / Box%size(:)
        
        exchange_energy = 0._DP
        
        do kz = 0, Box%wave(3)
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)
                    
                    wave_dot_orientation = dot_product(wave_vector, orientation_div_box)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                    
                    real_part = wave_dot_orientation * (wave_dot_orientation + 2._DP * &
                                real(this%structure(kx, ky, kz) * conjg(exp_Ikx), DP))

                    exchange_energy = exchange_energy + this%weight(kx, ky, kz) * real_part
                   
                end do
            
            end do
        
        end do
        
        exchange_energy = 4._DP*PI / product(Box%size) * exchange_energy

    end function Ewald_Summation_Reci_exchange_energy
    
    !> Exchange a particle -> update the ``structure factor''
    
    !> Add particle
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot+\vec{\mu}_{N+1}) e^{+i\vec{k}\cdot\vec{x}_{N+1}}
    !>  \f]
    !>
    
    !> Remove particle -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot-\vec{\mu}_{N}) e^{+i\vec{k}\cdot\vec{x}_{N}}
    !>  \f]
    !>
    
    pure subroutine Ewald_Summation_Reci_update_structure_exchange(this, Box, particle)

        class(Ewald_Summation_Reci), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange
        
        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: orientation_div_box
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx
        
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz
        
        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call fourier_i(Box%wave(3), position_div_box(3), exp_Ikx_3)
        
        orientation_div_box(:) = exchange_sign(particle%add) * particle%orientation(:) / Box%size(:)
        
        exchange = 0._DP
        
        do kz = 0, Box%wave(3)
            wave_vector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = real(kx, DP)
                    
                    wave_dot_orientation = dot_product(wave_vector, orientation_div_box)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                                                          
                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(wave_dot_orientation, 0._DP, DP) * exp_Ikx
                   
                end do
            
            end do
        
        end do
        
        exchange = 4._DP*PI / product(Box%size) * exchange

    end subroutine Ewald_Summation_Reci_update_structure_exchange

end module class_ewald_summation_reci
