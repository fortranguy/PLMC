module class_electronic_layer_correction

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use module_types_micro, only: Box_Parameters, Particle_Index
use module_physics_micro, only: fourier_i, Box_wave1_sym, exchange_sign, set_exp_kz
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Electronic_Layer_Correction

        private
        
        integer :: num_wave_vectors
        real(DP), dimension(:, :), allocatable :: wave_norm
        real(DP), dimension(:, :), allocatable :: weight
        complex(DP), dimension(:, :), allocatable :: structure_plus, structure_minus
        
    contains
    
        procedure :: construct => Electronic_Layer_Correction_construct
        procedure, private :: set_wave_norm => Electronic_Layer_Correction_set_wave_norm
        procedure, private :: set_weight => Electronic_Layer_Correction_set_weight
        procedure, private :: set_structure => Electronic_Layer_Correction_set_structure
        procedure :: destroy => Electronic_Layer_Correction_destroy
        procedure :: reset_structure => Electronic_Layer_Correction_reset_structure
        procedure, private :: get_structure_modulus => &
                              Electronic_Layer_Correction_get_structure_modulus
        procedure :: count_wave_vectors => Electronic_Layer_Correction_count_wave_vectors
        procedure :: total_energy => Electronic_Layer_Correction_total_energy
        procedure, private :: total_energy_structure => &
                              Electronic_Layer_Correction_total_energy_structure
        procedure, private :: total_energy_field => Electronic_Layer_Correction_total_energy_field
        
        procedure :: solo_field => Electronic_Layer_Correction_solo_field
        procedure :: test_field => Electronic_Layer_Correction_test_field
        
        procedure :: move_energy => Electronic_Layer_Correction_move_energy
        procedure :: update_structure_move => &
                     Electronic_Layer_Correction_update_structure_move
        procedure :: rotation_energy => Electronic_Layer_Correction_rotation_energy
        procedure :: update_structure_rotation => &
                     Electronic_Layer_Correction_update_structure_rotation
        procedure :: exchange_energy => Electronic_Layer_Correction_exchange_energy
        procedure :: update_structure_exchange => &
                     Electronic_Layer_Correction_update_structure_exchange
        
    end type Electronic_Layer_Correction
    
contains

    pure subroutine Electronic_Layer_Correction_construct(this, Box, this_spheres)
        
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        
        if (allocated(this%wave_norm)) deallocate(this%wave_norm)
        allocate(this%wave_norm(-Box%wave(1):Box%wave(1), &
                                -Box%wave(2):Box%wave(2)))
        
        if (allocated(this%weight)) deallocate(this%weight)
        allocate(this%weight(-Box%wave(1):Box%wave(1), &
                             -Box%wave(2):Box%wave(2)))
                                
        if (allocated(this%structure_plus)) deallocate(this%structure_plus)
        allocate(this%structure_plus(-Box%wave(1):Box%wave(1), &
                                     -Box%wave(2):Box%wave(2)))
        if (allocated(this%structure_minus)) deallocate(this%structure_minus)
        allocate(this%structure_minus(-Box%wave(1):Box%wave(1), &
                                      -Box%wave(2):Box%wave(2)))
                                      
        call this%set_wave_norm(Box)
        call this%set_weight(Box)
        call this%set_structure(Box, this_spheres)
    
    end subroutine Electronic_Layer_Correction_construct
    
    pure subroutine Electronic_Layer_Correction_set_wave_norm(this, Box)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
    
        integer :: kx, ky
        real(DP), dimension(2) :: wave_vector
        
        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
        
        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)
            
            this%wave_norm(kx, ky) = norm2(wave_vector)
            
        end do
        
        end do
    
    end subroutine Electronic_Layer_Correction_set_wave_norm
    
    !> \f[
    !>      w(\vec{k}^{2D}) = \frac{1}{k^{2D}(e^{k^{2D}L_z} - 1)}
    !> \f]
    
    pure subroutine Electronic_Layer_Correction_set_weight(this, Box)
        
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        
        integer :: kx, ky
        
        do ky = -Box%wave(2), Box%wave(2)
        
        do kx = -Box%wave(1), Box%wave(1)

            if (kx**2 + ky**2 /= 0) then
                this%weight(kx, ky) = 1._DP / this%wave_norm(kx, ky) / &
                                              (exp(this%wave_norm(kx, ky)*Box%size(3)) - 1._DP)
            else
                this%weight(kx, ky) = 0._DP
            end if

        end do
            
        end do
        
    end subroutine Electronic_Layer_Correction_set_weight
    
    !> Structure factor set :
    !> \f[
    !>      S_+(\vec{k}) = \sum_{i} (+\mu_{z,i}k^{2D} + i\vec{\mu}^{2D}_i\cdot\vec{k}^{2D})
    !>                     e^{+i\vec{k}^{2D}\cdot\vec{x}^{2D}_i} e^{+k^{2D}z_i}
    !> \f]
    !> We will also use a restricted definition later :
    !> \f[
    !>      S_{+, \underline{l}}(\vec{k}) = \sum_{i \neq l}
    !>                                      (+\mu_{z,i}k^{2D} + i\vec{\mu}^{2D}_i\cdot\vec{k}^{2D})
    !>                                      e^{+i\vec{k}^{2D}\cdot\vec{x}^{2D}_i} e^{+k^{2D}z_i}
    !> \f].

    pure subroutine Electronic_Layer_Correction_set_structure(this, Box, this_spheres)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres

        complex(DP) :: exp_Ikx
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kz_tab
        real(DP) :: exp_kz

        real(DP), dimension(2) :: position_div_box
        real(DP), dimension(2) :: wave_vector
        real(DP) :: wave_dot_orientation, wave_orientation_z
        integer :: kx, ky
        integer :: i_particle

        this%structure_plus(:, :) = cmplx(0._DP, 0._DP, DP)
        this%structure_minus(:, :) = cmplx(0._DP, 0._DP, DP)

        do i_particle = 1, this_spheres%get_num_particles()
        
            position_div_box(:) = 2._DP*PI * this_spheres%get_position_2d(i_particle)/Box%size(1:2)
            call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
            call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
            call set_exp_kz(Box%wave, this%wave_norm, this_spheres%get_position_z(i_particle), exp_kz_tab)
        
            do ky = -Box%wave(2), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
                
            do kx = -Box%wave(1), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)
            
                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_kz = exp_kz_tab(abs(kx), abs(ky))

                wave_dot_orientation = dot_product(wave_vector, &
                                                   this_spheres%get_orientation_2d(i_particle))
                wave_orientation_z = this%wave_norm(kx, ky) * this_spheres%get_orientation_z(i_particle)
                          
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                                              cmplx(+wave_orientation_z, wave_dot_orientation, DP) * &
                                              exp_Ikx * cmplx(exp_kz, 0._DP, DP)

                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                                               cmplx(-wave_orientation_z, wave_dot_orientation, DP) * &
                                               exp_Ikx / cmplx(exp_kz, 0._DP, DP)
            
            end do
            
            end do
            
        end do

    end subroutine Electronic_Layer_Correction_set_structure
    
    subroutine Electronic_Layer_Correction_destroy(this)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        
        if (allocated(this%wave_norm)) deallocate(this%wave_norm)
        if (allocated(this%weight)) deallocate(this%weight)
        if (allocated(this%structure_plus)) deallocate(this%structure_plus)
        if (allocated(this%structure_minus)) deallocate(this%structure_minus)
    
    end subroutine Electronic_Layer_Correction_destroy
    
    !> Reset the structure factor and print the drift
    
    subroutine Electronic_Layer_Correction_reset_structure(this, Box, this_spheres, i_step, modulus_unit)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: i_step
        integer, intent(in), optional :: modulus_unit

        real(DP), dimension(2) :: modulus_drifted, modulus_reset
        
        modulus_drifted(:) = this%get_structure_modulus(Box%wave)
        call this%set_structure(Box, this_spheres)
        modulus_reset(:) = this%get_structure_modulus(Box%wave)

        if (present(modulus_unit)) then
            write(modulus_unit, *) i_step, abs(modulus_reset(:) - modulus_drifted(:))
        end if
    
    end subroutine Electronic_Layer_Correction_reset_structure
    
    !> To calculate the drift of the structure factor

    pure function Electronic_Layer_Correction_get_structure_modulus(this, Box_wave) &
                  result(get_structure_modulus)

        class(Electronic_Layer_Correction), intent(in) :: this
        integer, dimension(:), intent(in) :: Box_wave
        real(DP), dimension(2) :: get_structure_modulus

        integer :: kx, ky

        get_structure_modulus(:) = 0._DP

        do ky = 0, Box_wave(2)
            do kx = -Box_wave1_sym(Box_wave, ky, 0), Box_wave(1)
                get_structure_modulus(1) = get_structure_modulus(1) + &
                                                     abs(this%structure_plus(kx, ky))
                get_structure_modulus(2) = get_structure_modulus(2) + &
                                                     abs(this%structure_minus(kx, ky))
            end do
        end do

    end function Electronic_Layer_Correction_get_structure_modulus
    
    ! Count the number of wave vectors

    subroutine Electronic_Layer_Correction_count_wave_vectors(this, Box_wave, wave_unit)

        class(Electronic_Layer_Correction), intent(inout) :: this
        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: wave_unit
        
        integer :: kx, ky

        this%num_wave_vectors = 0

        do ky = 0, Box_wave(2)
            do kx = -Box_wave1_sym(Box_wave, ky, 0), Box_wave(1)
                if (kx**2 + ky**2 /= 0) then

                    write(wave_unit, *) kx, ky
                    write(wave_unit, *)
                    write(wave_unit, *)

                    this%num_wave_vectors = this%num_wave_vectors + 1

                end if
            end do
        end do

    end subroutine Electronic_Layer_Correction_count_wave_vectors
    
    !> Total ELC energy
    
    pure function Electronic_Layer_Correction_total_energy(this, Box, using_field, this_spheres) &
                  result(total_energy)
    
        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        logical, intent(in), optional :: using_field
        type(Dipolar_Hard_Spheres), intent(in), optional :: this_spheres
        real(DP) :: total_energy
        
        if (present(using_field)) then
            if (using_field .and. present(this_spheres)) then
                total_energy = this%total_energy_field(Box, this_spheres)
            else
                total_energy = this%total_energy_structure(Box)
            end if
        else
            total_energy = this%total_energy_structure(Box)
        end if
    
    end function Electronic_Layer_Correction_total_energy
    
    pure function Electronic_Layer_Correction_total_energy_structure(this, Box) &
                  result(total_energy_structure)
        
        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP) :: total_energy_structure

        complex(DP) :: structure_product
        integer kx, ky

        total_energy_structure = 0._DP

        do ky = -Box%wave(2), Box%wave(2)
            do kx = -Box%wave(1), Box%wave(1)
                structure_product = this%structure_plus(kx, ky)*conjg(this%structure_minus(kx, ky))
                total_energy_structure = total_energy_structure + &
                                         this%weight(kx, ky) * 2._DP*real(structure_product, DP)
            end do
        end do
        
        total_energy_structure = PI / product(Box%size(1:2)) * total_energy_structure
        
    end function Electronic_Layer_Correction_total_energy_structure

    pure function Electronic_Layer_Correction_total_energy_field(this, Box, this_spheres) &
                  result(total_energy_field)
        
        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        real(DP) :: total_energy_field
        
        integer :: i_particle
        type(Particle_Index) :: particle

        total_energy_field = 0._DP

        do i_particle = 1, this_spheres%get_num_particles()
            particle%number = i_particle
            particle%position(:) = this_spheres%get_position(particle%number)
            particle%orientation(:) = this_spheres%get_orientation(particle%number)
            total_energy_field = total_energy_field + &
                                 dot_product(particle%orientation, &
                                             this%solo_field(Box, particle))
        end do
        
        total_energy_field = total_energy_field/2._DP
        
    end function Electronic_Layer_Correction_total_energy_field
    
    !> Field
    !> \f[
    !>      \vec{E}(\vec{r}_{N+1}) = \frac{2\pi}{S} \sum_{\vec{k}^{2D} \neq \vec{0}} w(k^{2D}) \{
    !>                                    \Re[ e^{(\vec{K}_+ \cdot \vec{r}_{N+1})} 
    !>                                              S_-^*(\vec{k}^{2D}) |\vec{K}_+) +
    !>                                         e^{(\vec{K}_-^* \cdot \vec{r}_{N+1})}
    !>                                              S_+(\vec{k}^{2D}) |\vec{K}_-^*) ] +
    !>                                    (\vec{\mu}_{N+1} \cdot \vec{K}_-^*)  |\vec{K}_+)
    !>                               \}  
    !> \f]
    
    pure function Electronic_Layer_Correction_solo_field(this, Box, particle) &
                  result(solo_field)
                  
        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP), dimension(num_dimensions) :: solo_field
        
        real(DP), dimension(2) :: position_div_box
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP) :: exp_Ikx
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kz_tab
        real(DP) :: exp_kz
        
        complex(DP), dimension(num_dimensions) :: wave_vector_plus, wave_vector_minus
        real(DP), dimension(2) :: wave_vector
        integer :: kx, ky
        
        position_div_box(:) = 2._DP*PI * particle%position(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_exp_kz(Box%wave, this%wave_norm, particle%position(3), exp_kz_tab)
        
        solo_field(:) = 0._DP
        
        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
        
        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)
            
            exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)
            exp_kz = exp_kz_tab(abs(kx), abs(ky))            
            
            wave_vector_plus(1:2) = cmplx(0._DP, wave_vector, DP)
            wave_vector_plus(3) = cmplx(+this%wave_norm(kx, ky), 0._DP, DP)
            
            wave_vector_minus(1:2) = cmplx(0._DP, wave_vector, DP)
            wave_vector_minus(3) = cmplx(-this%wave_norm(kx, ky), 0._DP, DP)
            
            solo_field(:) = solo_field(:) + this%weight(kx, ky) * &
                (real(conjg(this%structure_minus(kx, ky))*exp_Ikx*cmplx(exp_kz, 0._DP, DP) * &
                 wave_vector_plus(:) + &
                 this%structure_plus(kx, ky)*conjg(exp_Ikx)/cmplx(exp_kz, 0._DP, DP) * &
                 conjg(wave_vector_minus(:)) + &
                 dot_product(particle%orientation, conjg(wave_vector_minus)) * &
                 wave_vector_plus(:), DP))
            
        end do
        
        end do
        
        solo_field(:) = 2._DP*PI / product(Box%size(1:2)) * solo_field(:)        
                  
    end function Electronic_Layer_Correction_solo_field

    !> Field
    !> \f[
    !>      \vec{E}(\vec{r}_{N+1}) =-\frac{2\pi}{S} \sum_{\vec{k}^{2D} \neq \vec{0}} w(k^{2D}) \{
    !>                                    \Re[ e^{(\vec{K}_+ \cdot \vec{r}_{N+1})} 
    !>                                              S_-^*(\vec{k}^{2D}) |\vec{K}_+) +
    !>                                         e^{(\vec{K}_-^* \cdot \vec{r}_{N+1})}
    !>                                              S_+(\vec{k}^{2D}) |\vec{K}_-^*) ] +
    !>                                    (\vec{\mu}_{N+1} \cdot \vec{K}_-^*)  |\vec{K}_+)
    !>                               \}  
    !> \f]

    pure function Electronic_Layer_Correction_test_field(this, Box, particle) &
                  result(test_field)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP), dimension(num_dimensions) :: test_field

        test_field(:) = -this%solo_field(Box, particle)

    end function Electronic_Layer_Correction_test_field
    
    !> Move

    !> Difference of Energy
    !>  \f[ \Delta U = \frac{\pi}{S} \sum_{\vec{k}^{2D} \neq \vec{0}} w(\vec{k}^{2D}) \Delta S^2
    !>  \f]
    !> \f[
    !>  \Delta S^2 = 2\Re[
    !>                   S_+^\prime S_-^{\prime*} - S_+ S_-^*
    !>               ]
    !> \f]
    
    !> Algebra & Notation
    !> \f[
    !>  \Delta S^2 = 2\Re(
    !>                   S_{+,\underline{l}} (s_{-,l}^{\prime*} - s_{-,l}^*) +
    !>                   S_{-,\underline{l}}^* (s_{+,l}^{\prime} - s_{+,l})
    !>               )
    !> \f]

    pure function Electronic_Layer_Correction_move_energy(this, Box, old, new) result(move_energy)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: move_energy

        real(DP), dimension(2) :: new_position_div_box, old_position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP) :: exp_IkxNew
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzNew_tab
        real(DP) :: exp_kzNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP) :: exp_IkxOld
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzOld_tab
        real(DP) :: exp_kzOld

        complex(DP) :: structure_i, delta_structure_i
        
        real(DP) :: real_part1, real_part2
        
        real(DP), dimension(2) :: wave_vector
        real(DP) :: wave_dot_orientation, wave_orientation_z
        integer :: kx, ky

        new_position_div_box(:) = 2._DP*PI * new%position(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call set_exp_kz(Box%wave, this%wave_norm, new%position(3), exp_kzNew_tab)
        
        old_position_div_box(:) = 2._DP*PI * old%position(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call set_exp_kz(Box%wave, this%wave_norm, old%position(3), exp_kzOld_tab)

        move_energy = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)
                
                wave_orientation_z = this%wave_norm(kx, ky) * new%orientation(3)
                wave_dot_orientation = dot_product(wave_vector, new%orientation(1:2))
                
                exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky)
                exp_kzNew = exp_kzNew_tab(abs(kx), ky)
                exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky)
                exp_kzOld = exp_kzOld_tab(abs(kx), ky)

                ! S_{+,\underline{l}} (s_{-,l}^{\prime*} - s_{-,l}^*) +

                structure_i = cmplx(+wave_orientation_z, wave_dot_orientation, DP) * &
                              exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP)

                delta_structure_i = cmplx(-wave_orientation_z, wave_dot_orientation, DP) * &
                (exp_IkxNew/cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP))

                real_part1 = real((this%structure_plus(kx, ky) - structure_i) * &
                            conjg(delta_structure_i), DP)

                ! S_{-,\underline{l}}^* (s_{+,l}^{\prime} - s_{+,l})
                
                structure_i = cmplx(-wave_orientation_z, wave_dot_orientation, DP) * &
                              exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP)

                delta_structure_i = cmplx(+wave_orientation_z, wave_dot_orientation, DP) * &
                (exp_IkxNew*cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP))

                real_part2 = real(conjg(this%structure_minus(kx, ky) - structure_i) * &
                            delta_structure_i, DP)
                            
                !   Accumulation
                
                move_energy = move_energy + this%weight(kx, ky) * 2._DP * (real_part1 + real_part2)

            end do
        
        end do

        move_energy = 2._DP*PI / product(Box%size(1:2)) * move_energy

    end function Electronic_Layer_Correction_move_energy

    !> Update position -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = (\pm k^{2D}\mu_{l ,z} + i(\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))
    !>                       (e^{\pm k^{2D}z_l^\prime} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D\prime}_l)} -
    !>                        e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)})
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structure_move(this, Box, old, new)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        
        real(DP), dimension(2) :: new_position_div_box, old_position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP) :: exp_IkxNew
        
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzNew_tab
        real(DP) :: exp_kzNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP) :: exp_IkxOld
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzOld_tab
        real(DP) :: exp_kzOld

        real(DP), dimension(2) :: wave_vector
        real(DP) :: wave_dot_orientation, wave_orientation_z
        integer :: kx, ky

        new_position_div_box(:) = 2._DP*PI * new%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call set_exp_kz(Box%wave, this%wave_norm, new%position(3), exp_kzNew_tab)
        
        old_position_div_box(:) = 2._DP*PI * old%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call set_exp_kz(Box%wave, this%wave_norm, old%position(3), exp_kzOld_tab)

        do ky = 0, Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                wave_dot_orientation = dot_product(wave_vector, new%orientation(1:2))
                wave_orientation_z = this%wave_norm(kx, ky) * new%orientation(3)

                exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky)
                exp_kzNew = exp_kzNew_tab(abs(kx), ky)
                
                exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky)
                exp_kzOld = exp_kzOld_tab(abs(kx), ky)
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                    cmplx(+wave_orientation_z, wave_dot_orientation, DP) * &
                    (exp_IkxNew*cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP))
                                              
                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                    cmplx(-wave_orientation_z, wave_dot_orientation, DP) * &
                    (exp_IkxNew/cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP))

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structure_move
    
    !> Rotate

    !> Difference of Energy
    !>  \f[ \Delta U = \frac{\pi}{S} \sum_{\vec{k}^{2D} \neq \vec{0}} w(\vec{k}^{2D}) \Delta S^2
    !>  \f]
    !> \f[
    !>  \Delta S^2 = 2\Re[
    !>                   S_+^\prime S_-^{\prime*} - S_+ S_-^*
    !>               ]
    !> \f]
    
    !> Algebra & Notation
    !> \f[
    !>  \Delta S^2 = 2\Re(
    !>                   s_{+, l}^{\prime} s_{-, l}^{\prime*} - s_{+, l} s_{-, l}^* +
    !>                   S_{+, \underline{l}} (s_{-, l}^{\prime*} - s_{-, l}^*) +
    !>                   S_{-, \underline{l}}^* (s_{+, l}^{\prime} - s_{+, l})
    !>               )
    !> \f]
    
    pure function Electronic_Layer_Correction_rotation_energy(this, Box, old, new) &
         result(rotation_energy)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation_energy

        real(DP), dimension(2) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP) :: exp_Ikx
        
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_Ikz_tab
        real(DP) :: exp_Ikz

        complex(DP) :: structure_i, delta_structure_i
        real(DP) :: real_part0, real_part1, real_part2
        
        real(DP), dimension(2) :: wave_vector
        real(DP) :: kMnew_z, kMold_z
        real(DP) :: k_dot_mNew, k_dot_mOld
        integer :: kx, ky

        position_div_box(:) = 2._DP*PI * new%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_exp_kz(Box%wave, this%wave_norm, new%position(3), exp_Ikz_tab)

        rotation_energy = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                kMnew_z = this%wave_norm(kx, ky) * new%orientation(3)
                k_dot_mNew = dot_product(wave_vector, new%orientation(1:2))
                
                kMold_z = this%wave_norm(kx, ky) * old%orientation(3)
                k_dot_mOld = dot_product(wave_vector, old%orientation(1:2))

                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_Ikz = exp_Ikz_tab(abs(kx), ky)
                
                ! s_{+,l}^{\prime} s_{-,l}^{\prime*} - s_{+,l} s_{-,l}^*
                
                real_part0 = (-kMnew_z**2 + k_dot_mNew**2) + (kMold_z**2 - k_dot_mOld**2)

                ! S_{+,\underline{l}} (s_{-,l}^{\prime*} - s_{-,l}^*)

                structure_i = cmplx(+kMold_z, k_dot_mOld, DP) * exp_Ikx*cmplx(exp_Ikz, 0._DP, DP)

                delta_structure_i = cmplx(-kMnew_z, k_dot_mNew, DP) - cmplx(-kMold_z, k_dot_mOld, DP)
                delta_structure_i = delta_structure_i * exp_Ikx / cmplx(exp_Ikz, 0._DP, DP)

                real_part1 = real((this%structure_plus(kx, ky) - structure_i) * &
                            conjg(delta_structure_i), DP)

                ! S_{-,\underline{l}}^* (s_{+,l}^{\prime} - s_{+,l})
                                 
                structure_i = cmplx(-kMold_z, k_dot_mOld, DP) * exp_Ikx / cmplx(exp_Ikz, 0._DP, DP)

                delta_structure_i = cmplx(+kMnew_z, k_dot_mNew, DP) - cmplx(+kMold_z, k_dot_mOld, DP)
                delta_structure_i = delta_structure_i * exp_Ikx * cmplx(exp_Ikz, 0._DP, DP)

                real_part2 = real(conjg(this%structure_minus(kx, ky) - structure_i) * &
                            delta_structure_i, DP)
                            
                !   Accumulation

                rotation_energy = rotation_energy + this%weight(kx, ky) * 2._DP * &
                                  (real_part0 + real_part1 + real_part2)

            end do

        end do

        rotation_energy = 2._DP*PI / product(Box%size(1:2)) * rotation_energy

    end function Electronic_Layer_Correction_rotation_energy

    !> Update moment -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = [\pm (k^{2D}\mu_{l ,z}^\prime - k^{2D}\mu_{l ,z}) +
    !>                       i((\vec{k}^{2D}\cdot\vec{\mu}^{2D\prime}_l) -
    !>                         (\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))]
    !>                       e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)}
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structure_rotation(this, Box, old, new)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new

        real(DP), dimension(2) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP) :: exp_Ikx
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_Ikz_tab
        real(DP) :: exp_Ikz

        real(DP), dimension(2) :: wave_vector
        real(DP) :: kdeltaMcol_z, k_dot_deltaMcol
        integer :: kx, ky

        position_div_box(:) = 2._DP*PI * new%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_exp_kz(Box%wave, this%wave_norm, new%position(3), exp_Ikz_tab)

        do ky = 0, Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                exp_Ikz = exp_Ikz_tab(abs(kx), ky)
                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)

                kdeltaMcol_z = this%wave_norm(kx, ky) * (new%orientation(3) - old%orientation(3))
                k_dot_deltaMcol = dot_product(wave_vector, &
                                              new%orientation(1:2) - old%orientation(1:2))
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                    cmplx(+kdeltaMcol_z, k_dot_deltaMcol, DP) * &
                    exp_Ikx * cmplx(exp_Ikz, 0._DP, DP)

                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                    cmplx(-kdeltaMcol_z, k_dot_deltaMcol, DP) * &
                    exp_Ikx / cmplx(exp_Ikz, 0._DP, DP)

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structure_rotation
    
    !> Energy of 1 dipole with others
    
    !> Addition :
    
    !> Difference of Energy
    !>  \f[ \Delta U = \frac{\pi}{S} \sum_{\vec{k}^{2D} \neq \vec{0}} w(\vec{k}^{2D}) \Delta S^2
    !>  \f]
    !> \f[
    !>  \Delta S^2 = 2\Re[S_+^{N+1} S_-^{N+1*} - S_+^N S_-^{N*}]
    !> \f]
    
    !> Algebra & Notation
    !> \f[
    !>  \Delta S^2 = 2\Re(
    !>                   s_{+, N+1} s_{-, N+1}^* +
    !>                   S_+^N s_{-, N+1}^* + S_{-}^{N*} s_{+, N+1}
    !>               )
    !> \f]
    
    !> Summary: only the sign of \f[\vec{\mu}\f] changes.

    pure function Electronic_Layer_Correction_exchange_energy(this, Box, particle) result(exchange_energy)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange_energy
        
        real(DP), dimension(2) :: position_div_box
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP) :: exp_Ikx
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_Ikz_tab
        real(DP) :: exp_Ikz
        
        complex(DP) :: structure_i
        real(DP) :: real_part0, real_part1, real_part2
        
        real(DP), dimension(2) :: wave_vector
        real(DP) :: wave_orientation_z, wave_dot_orientation
        integer :: kx, ky
        real(DP) :: exchg_sign
        
        position_div_box(:) = 2._DP*PI * particle%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_exp_kz(Box%wave, this%wave_norm, particle%position(3), exp_Ikz_tab)
        
        exchg_sign = exchange_sign(particle%add)
        
        exchange_energy = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)
                
                wave_orientation_z = this%wave_norm(kx, ky) * exchg_sign * particle%orientation(3)
                wave_dot_orientation = dot_product(wave_vector, &
                                                   exchg_sign * particle%orientation(1:2))
                                                
                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_Ikz = exp_Ikz_tab(abs(kx), ky)
                
                ! s_{+, N+1} s_{-, N+1}^*
                
                real_part0 = -wave_orientation_z**2 + wave_dot_orientation**2
                
                ! S_+^N s_{-, N+1}^*
                
                structure_i = cmplx(-wave_orientation_z, wave_dot_orientation, DP) * exp_Ikx/cmplx(exp_Ikz, 0._DP, DP)
                real_part1 = real(this%structure_plus(kx, ky) * conjg(structure_i), DP)
                
                ! S_{-}^{N*} s_{+, N+1}
                
                structure_i = cmplx(wave_orientation_z, wave_dot_orientation, DP) * exp_Ikx*cmplx(exp_Ikz, 0._DP, DP)
                real_part2 = real(conjg(this%structure_minus(kx, ky)) * structure_i, DP)
                
                !   Accumulation

                exchange_energy = exchange_energy + this%weight(kx, ky) * &
                                                    2._DP * (real_part0 + real_part1 + real_part2)
                
            end do
        
        end do

        exchange_energy = 2._DP * PI/product(Box%size(1:2)) * exchange_energy

    end function Electronic_Layer_Correction_exchange_energy
    
    !> Exchange a particle -> update the ``structure factor''
    
    !> Add particle
    !>  \f[
    !>      \Delta S_{\pm} = [\pm k^{2D}\mu_{N+1 ,z} + i(\vec{k}^{2D}\cdot\vec{\mu}^{2D}_{N+1})]
    !>                       e^{\pm k^{2D}z_{N+1}} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_{N+1})}
    !>  \f]
    !>

    pure subroutine Electronic_Layer_Correction_update_structure_exchange(this, Box, particle)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        
        real(DP), dimension(2) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP) :: exp_Ikx
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_Ikz_tab
        real(DP) :: exp_Ikz

        real(DP), dimension(2) :: wave_vector
        real(DP) :: wave_orientation_z, wave_dot_orientation
        integer :: kx, ky
        real(DP) :: exchg_sign

        position_div_box(:) = 2._DP*PI * particle%position(1:2) / Box%size(1:2)
        call fourier_i(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call fourier_i(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_exp_kz(Box%wave, this%wave_norm, particle%position(3), exp_Ikz_tab)
        
        exchg_sign = exchange_sign(particle%add)

        do ky = 0, Box%wave(2)

            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)

                wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                wave_orientation_z = this%wave_norm(kx, ky) * exchg_sign * particle%orientation(3)
                wave_dot_orientation = dot_product(wave_vector, &
                                                   exchg_sign * particle%orientation(1:2))
                
                exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_Ikz = exp_Ikz_tab(abs(kx), ky)
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                                              cmplx(+wave_orientation_z, wave_dot_orientation, DP) * &
                                              exp_Ikx * cmplx(exp_Ikz, 0._DP, DP)
                
                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                                               cmplx(-wave_orientation_z, wave_dot_orientation, DP) * &
                                               exp_Ikx / cmplx(exp_Ikz, 0._DP, DP)

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structure_exchange

end module class_electronic_layer_correction
