module class_dipolarSpheres

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
        procedure, private :: set_structures => Electronic_Layer_Correction_set_structures
        procedure :: destroy => Electronic_Layer_Correction_destroy
        procedure :: reset_structures => Electronic_Layer_Correction_reset_structures
        procedure, private :: get_structures_modulus => &
                              Electronic_Layer_Correction_get_structures_modulus
        procedure :: count_wave_vectors => Electronic_Layer_Correction_count_wave_vectors
        procedure, private :: total => Electronic_Layer_Correction_total
        
        !>     ELC: delta
        procedure :: move => Electronic_Layer_Correction_move
        procedure :: update_structures_move => &
                     Electronic_Layer_Correction_update_structures_move
        procedure :: rotation => Electronic_Layer_Correction_rotation
        procedure :: update_structures_rotate => &
                     Electronic_Layer_Correction_update_structures_rotate
        procedure :: exchange => Electronic_Layer_Correction_exchange
        procedure :: update_structures_exchange => &
                     Electronic_Layer_Correction_update_structures_exchange
        !>     ELC: total
        
        
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
        call this%set_structures(Box, this_spheres)
    
    end subroutine Electronic_Layer_Correction_construct
    
    pure subroutine Electronic_Layer_Correction_set_wave_norm(this, Box)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
    
        integer :: kx, ky
        real(DP), dimension(num_dimensions-1) :: wave_vector
        
        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = real(ky, DP)
        
        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = real(kx, DP)
            
            this%wave_norm(kx, ky) = 2._DP*PI * norm2(wave_vector(:)/Box%size(1:num_dimensions-1))
            
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
    !> We will also use a restricted defsetion later :
    !> \f[
    !>      S_{+, \underline{l}}(\vec{k}) = \sum_{i \neq l}
    !>                                      (+\mu_{z,i}k^{2D} + i\vec{\mu}^{2D}_i\cdot\vec{k}^{2D})
    !>                                      e^{+i\vec{k}^{2D}\cdot\vec{x}^{2D}_i} e^{+k^{2D}z_i}
    !> \f].

    pure subroutine Electronic_Layer_Correction_set_structures(this, Box, this_spheres)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        real(DP) :: exp_kzCol

        real(DP), dimension(num_dimensions-1) :: xColOverL
        real(DP), dimension(num_dimensions-1) :: mColOverL
        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky
        integer :: iCol

        this%structure_plus(:, :) = cmplx(0._DP, 0._DP, DP)
        this%structure_minus(:, :) = cmplx(0._DP, 0._DP, DP)

        do iCol = 1, this_spheres%get_num_particles()
        
            xColOverL(:) = 2._DP*PI * this_spheres%get_position_2d(iCol)/Box%size(1:2)
            call fourier_i(Box%wave(1), xColOverL(1), exp_Ikx_1)
            call fourier_i(Box%wave(2), xColOverL(2), exp_Ikx_2)
            
            mColOverL(:) = 2._DP*PI * this_spheres%get_orientation_2d(iCol)/Box%size(1:2)
        
            do ky = -Box%wave(2), Box%wave(2)

                wave_vector(2) = real(ky, DP)

            do kx = -Box%wave(1), Box%wave(1)

                wave_vector(1) = real(kx, DP)
            
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_kzCol = exp(this%wave_norm(kx, ky) * this_spheres%get_position_z(iCol))

                k_dot_mCol = dot_product(wave_vector, mColOverL)
                kMcol_z = this%wave_norm(kx, ky) * this_spheres%get_orientation_z(iCol)
                          
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                                              cmplx(+kMcol_z, k_dot_mCol, DP) * &
                                              exp_IkxCol * cmplx(exp_kzCol, 0._DP, DP)

                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                                               cmplx(-kMcol_z, k_dot_mCol, DP) * &
                                               exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)
            
            end do
            
            end do
            
        end do

    end subroutine Electronic_Layer_Correction_set_structures
    
    subroutine Electronic_Layer_Correction_destroy(this)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        
        if (allocated(this%wave_norm)) deallocate(this%wave_norm)
        if (allocated(this%weight)) deallocate(this%weight)
        if (allocated(this%structure_plus)) deallocate(this%structure_plus)
        if (allocated(this%structure_minus)) deallocate(this%structure_minus)
    
    end subroutine Electronic_Layer_Correction_destroy
    
    !> Reset the structure factor and print the drift
    
    subroutine Electronic_Layer_Correction_reset_structures(this, Box, this_spheres, iStep, modulus_unit)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Dipolar_Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP), dimension(2) :: modulus_drifted, modulus_reset
        
        modulus_drifted(:) = this%get_structures_modulus(Box%wave)
        call this%set_structures(Box, this_spheres)
        modulus_reset(:) = this%get_structures_modulus(Box%wave)
        
        write(modulus_unit, *) iStep, abs(modulus_reset(:) - modulus_drifted(:))
    
    end subroutine Electronic_Layer_Correction_reset_structures
    
    !> To calculate the drift of the structure factor

    pure function Electronic_Layer_Correction_get_structures_modulus(this, Box_wave) &
                  result(get_structures_modulus)

        class(Electronic_Layer_Correction), intent(in) :: this
        integer, dimension(:), intent(in) :: Box_wave
        real(DP), dimension(2) :: get_structures_modulus

        integer :: kx, ky

        get_structures_modulus(:) = 0._DP

        do ky = 0, Box_wave(2)
            do kx = -Box_wave1_sym(Box_wave, ky, 0), Box_wave(1)
                get_structures_modulus(1) = get_structures_modulus(1) + &
                                                     abs(this%structure_plus(kx, ky))
                get_structures_modulus(2) = get_structures_modulus(2) + &
                                                     abs(this%structure_minus(kx, ky))
            end do
        end do

    end function Electronic_Layer_Correction_get_structures_modulus
    
    ! Count the number of wave vectors

    subroutine Electronic_Layer_Correction_count_wave_vectors(this, Box_wave, wave_vectors_unit)

        class(Electronic_Layer_Correction), intent(inout) :: this
        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: wave_vectors_unit
        
        integer :: kx, ky

        this%num_wave_vectors = 0

        do ky = 0, Box_wave(2)
            do kx = -Box_wave1_sym(Box_wave, ky, 0), Box_wave(1)
                if (kx**2 + ky**2 /= 0) then
                    write(wave_vectors_unit, *) kx, ky
                    write(wave_vectors_unit, *)
                    write(wave_vectors_unit, *)

                    this%num_wave_vectors = this%num_wave_vectors + 1

                end if
            end do
        end do

    end subroutine Electronic_Layer_Correction_count_wave_vectors
    
    !> Total ELC energy
    
    pure function Electronic_Layer_Correction_total(this, Box) result(total)
        
        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP) :: total

        complex(DP) :: structure_product
        integer kx, ky

        total = 0._DP

        do ky = -Box%wave(2), Box%wave(2)
            do kx = -Box%wave(1), Box%wave(1)
                structure_product = this%structure_plus(kx, ky)*conjg(this%structure_minus(kx, ky))
                total = total + this%weight(kx, ky) * 2._DP*real(structure_product, DP)
            end do
        end do
        
        total = PI/product(Box%size(1:2)) * total
        
    end function Electronic_Layer_Correction_total
    
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

    pure function Electronic_Layer_Correction_move(this, Box, xOld, xNew, mCol) result(move)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: move

        real(DP), dimension(num_dimensions-1) :: xNewOverL, xOldOverL
        real(DP), dimension(num_dimensions-1) :: mColOverL

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
        
        real(DP) :: realPart1, realPart2
        
        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky

        xNewOverL(:) = 2._DP*PI * xNew(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        call fourier_i(Box%wave(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), xNewOverL(2), exp_IkxNew_2)
        call set_exp_kz(Box%wave, this%wave_norm, xNew(3), exp_kzNew_tab)
        
        xOldOverL(:) = 2._DP*PI * xOld(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        call fourier_i(Box%wave(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), xOldOverL(2), exp_IkxOld_2)
        call set_exp_kz(Box%wave, this%wave_norm, xOld(3), exp_kzOld_tab)

        mColOverL(:) = 2._DP*PI * mCol(1:num_dimensions-1)/Box%size(1:num_dimensions-1)

        move = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = real(ky, DP)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = real(kx, DP)
                
                kMcol_z = this%wave_norm(kx, ky) * mCol(3)
                k_dot_mCol = dot_product(wave_vector, mColOverL)
                
                exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky)
                exp_kzNew = exp_kzNew_tab(abs(kx), ky)
                exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky)
                exp_kzOld = exp_kzOld_tab(abs(kx), ky)

                ! S_{+,\underline{l}} (s_{-,l}^{\prime*} - s_{-,l}^*) +

                structure_i = cmplx(+kMcol_z, k_dot_mCol, DP) * exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP)

                delta_structure_i = cmplx(-kMcol_z, k_dot_mCol, DP) * &
                (exp_IkxNew/cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP))

                realPart1 = real((this%structure_plus(kx, ky) - structure_i) * &
                            conjg(delta_structure_i), DP)

                ! S_{-,\underline{l}}^* (s_{+,l}^{\prime} - s_{+,l})
                
                structure_i = cmplx(-kMcol_z, k_dot_mCol, DP) * exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP)

                delta_structure_i = cmplx(+kMcol_z, k_dot_mCol, DP) * &
                (exp_IkxNew*cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP))

                realPart2 = real(conjg(this%structure_minus(kx, ky) - structure_i) * &
                            delta_structure_i, DP)
                            
                !   Accumulation
                
                move = move + this%weight(kx, ky) * 2._DP * &
                                     (realPart1 + realPart2)

            end do
        
        end do

        move = 2._DP * PI/product(Box%size(1:2)) * move

    end function Electronic_Layer_Correction_move

    !> Update position -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = (\pm k^{2D}\mu_{l ,z} + i(\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))
    !>                       (e^{\pm k^{2D}z_l^\prime} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D\prime}_l)} -
    !>                        e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)})
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structures_move(this, Box, xOld, xNew, mCol)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        
        real(DP), dimension(num_dimensions-1) :: xNewOverL, xOldOverL
        real(DP), dimension(num_dimensions-1) :: mColOverL

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

        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky

        xNewOverL(:) = 2._DP*PI * xNew(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        call fourier_i(Box%wave(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), xNewOverL(2), exp_IkxNew_2)
        call set_exp_kz(Box%wave, this%wave_norm, xNew(3), exp_kzNew_tab)
        
        xOldOverL(:) = 2._DP*PI * xOld(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        call fourier_i(Box%wave(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), xOldOverL(2), exp_IkxOld_2)
        call set_exp_kz(Box%wave, this%wave_norm, xOld(3), exp_kzOld_tab)

        mColOverL(:) = 2._DP*PI * mCol(1:num_dimensions-1)/Box%size(1:num_dimensions-1)

        do ky = 0, Box%wave(2)
            wave_vector(2) = real(ky, DP)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = real(kx, DP)

                k_dot_mCol = dot_product(wave_vector, mColOverL)
                kMcol_z = this%wave_norm(kx, ky) * mCol(3)

                exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky)
                exp_kzNew = exp_kzNew_tab(abs(kx), ky)
                
                exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky)
                exp_kzOld = exp_kzOld_tab(abs(kx), ky)
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                    cmplx(+kMcol_z, k_dot_mCol, DP) * &
                    (exp_IkxNew*cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld*cmplx(exp_kzOld, 0._DP, DP))
                                              
                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                    cmplx(-kMcol_z, k_dot_mCol, DP) * &
                    (exp_IkxNew/cmplx(exp_kzNew, 0._DP, DP) - exp_IkxOld/cmplx(exp_kzOld, 0._DP, DP))

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structures_move
    
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
    
    pure function Electronic_Layer_Correction_rotation(this, Box, xCol, mOld, mNew) result(rotation)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: rotation

        real(DP), dimension(num_dimensions-1) :: xColOverL
        real(DP), dimension(num_dimensions-1) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol

        complex(DP) :: structure_i, delta_structure_i
        real(DP) :: realPart0, realPart1, realPart2
        
        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: kMnew_z, kMold_z
        real(DP) :: k_dot_mNew, k_dot_mOld
        integer :: kx, ky

        xColOverL(:) = 2._DP*PI * xCol(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(Box%wave, this%wave_norm, xCol(3), exp_kzCol_tab)

        mNewOverL(:) = 2._DP*PI * mNew(1:2)/Box%size(1:2)
        mOldOverL(:) = 2._DP*PI * mOld(1:2)/Box%size(1:2)

        rotation = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = real(ky, DP)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = real(kx, DP)

                kMnew_z = this%wave_norm(kx, ky) * mNew(3)
                k_dot_mNew = dot_product(wave_vector, mNewOverL)
                
                kMold_z = this%wave_norm(kx, ky) * mOld(3)
                k_dot_mOld = dot_product(wave_vector, mOldOverL)

                exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky)
                exp_kzCol = exp_kzCol_tab(abs(kx), ky)
                
                ! s_{+,l}^{\prime} s_{-,l}^{\prime*} - s_{+,l} s_{-,l}^*
                
                realPart0 = (-kMnew_z**2 + k_dot_mNew**2) + (kMold_z**2 - k_dot_mOld**2)

                ! S_{+,\underline{l}} (s_{-,l}^{\prime*} - s_{-,l}^*)

                structure_i = cmplx(+kMold_z, k_dot_mOld, DP) * exp_IkxCol*cmplx(exp_kzCol, 0._DP, DP)

                delta_structure_i = cmplx(-kMnew_z, k_dot_mNew, DP) - cmplx(-kMold_z, k_dot_mOld, DP)
                delta_structure_i = delta_structure_i * exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)

                realPart1 = real((this%structure_plus(kx, ky) - structure_i) * &
                            conjg(delta_structure_i), DP)

                ! S_{-,\underline{l}}^* (s_{+,l}^{\prime} - s_{+,l})
                                 
                structure_i = cmplx(-kMold_z, k_dot_mOld, DP) * exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)

                delta_structure_i = cmplx(+kMnew_z, k_dot_mNew, DP) - cmplx(+kMold_z, k_dot_mOld, DP)
                delta_structure_i = delta_structure_i * exp_IkxCol * cmplx(exp_kzCol, 0._DP, DP)

                realPart2 = real(conjg(this%structure_minus(kx, ky) - structure_i) * &
                            delta_structure_i, DP)
                            
                !   Accumulation

                rotation = rotation + this%weight(kx, ky) * 2._DP * &
                                     (realPart0 + realPart1 + realPart2)

            end do

        end do

        rotation = 2._DP * PI/product(Box%size(1:2)) * rotation

    end function Electronic_Layer_Correction_rotation

    !> Update moment -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = [\pm (k^{2D}\mu_{l ,z}^\prime - k^{2D}\mu_{l ,z}) +
    !>                       i((\vec{k}^{2D}\cdot\vec{\mu}^{2D\prime}_l) -
    !>                         (\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))]
    !>                       e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)}
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structures_rotate(this, Box, xCol, mOld, mNew)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew

        real(DP), dimension(num_dimensions-1) :: xColOverL
        real(DP), dimension(num_dimensions-1) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol

        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: kdeltaMcol_z, k_dot_deltaMcol
        integer :: kx, ky

        xColOverL(:) = 2._DP*PI * xCol(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(Box%wave, this%wave_norm, xCol(3), exp_kzCol_tab)

        mNewOverL(:) = 2._DP*PI * mNew(1:num_dimensions-1)/Box%size(1:num_dimensions-1)
        mOldOverL(:) = 2._DP*PI * mOld(1:num_dimensions-1)/Box%size(1:num_dimensions-1)

        do ky = 0, Box%wave(2)
            wave_vector(2) = real(ky, DP)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = real(kx, DP)

                exp_kzCol = exp_kzCol_tab(abs(kx), ky)
                exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky)

                kdeltaMcol_z = this%wave_norm(kx, ky) * (mNew(3) - mOld(3))
                k_dot_deltaMcol = dot_product(wave_vector, mNewOverL - mOldOverL)
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                    cmplx(+kdeltaMcol_z, k_dot_deltaMcol, DP) * &
                    exp_IkxCol * cmplx(exp_kzCol, 0._DP, DP)

                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                    cmplx(-kdeltaMcol_z, k_dot_deltaMcol, DP) * &
                    exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structures_rotate
    
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

    pure function Electronic_Layer_Correction_exchange(this, Box, xCol, mCol) result(exchange)

        class(Electronic_Layer_Correction), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: exchange
        
        real(DP), dimension(num_dimensions-1) :: xColOverL
        real(DP), dimension(num_dimensions-1) :: mColOverL
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol
        
        complex(DP) :: structure_i
        real(DP) :: realPart0, realPart1, realPart2
        
        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: kMcol_z, k_dot_mCol
        integer :: kx, ky
        
        xColOverL(:) = 2._DP*PI * xCol(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(Box%wave, this%wave_norm, xCol(3), exp_kzCol_tab)
        
        mColOverL(:) = 2._DP*PI * mCol(1:2)/Box%size(1:2)
        
        exchange = 0._DP

        do ky = 0, Box%wave(2)
            wave_vector(2) = real(ky, DP)
        
            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)
                wave_vector(1) = real(kx, DP)
                
                kMcol_z = this%wave_norm(kx, ky) * mCol(3)
                k_dot_mCol = dot_product(wave_vector, mColOverL)
                                                
                exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky)
                exp_kzCol = exp_kzCol_tab(abs(kx), ky)
                
                ! s_{+, N+1} s_{-, N+1}^*
                
                realPart0 = -kMcol_z**2 + k_dot_mCol**2
                
                ! S_+^N s_{-, N+1}^*
                
                structure_i = cmplx(-kMcol_z, k_dot_mCol, DP) * exp_IkxCol/cmplx(exp_kzCol, 0._DP, DP)
                realPart1 = real(this%structure_plus(kx, ky) * conjg(structure_i), DP)
                
                ! S_{-}^{N*} s_{+, N+1}
                
                structure_i = cmplx(kMcol_z, k_dot_mCol, DP) * exp_IkxCol*cmplx(exp_kzCol, 0._DP, DP)
                realPart2 = real(conjg(this%structure_minus(kx, ky)) * structure_i, DP)
                
                !   Accumulation

                exchange = exchange + this%weight(kx, ky) * &
                                         2._DP * (realPart0 + realPart1 + realPart2)
                
            end do
        
        end do

        exchange = 2._DP * PI/product(Box%size(1:2)) * exchange

    end function Electronic_Layer_Correction_exchange
    
    !> Exchange a particle -> update the ``structure factor''
    
    !> Add particle
    !>  \f[
    !>      \Delta S_{\pm} = [\pm k^{2D}\mu_{N+1 ,z} + i(\vec{k}^{2D}\cdot\vec{\mu}^{2D}_{N+1})]
    !>                       e^{\pm k^{2D}z_{N+1}} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_{N+1})}
    !>  \f]
    !>

    pure subroutine Electronic_Layer_Correction_update_structures_exchange(this, Box, xCol, mCol)

        class(Electronic_Layer_Correction), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mCol
        
        real(DP), dimension(num_dimensions-1) :: xColOverL
        real(DP), dimension(num_dimensions-1) :: mColOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Box%wave(1), 0:Box%wave(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol

        real(DP), dimension(num_dimensions-1) :: wave_vector
        real(DP) :: kMcol_z, k_dot_mCol
        integer :: kx, ky

        xColOverL(:) = 2._DP*PI * xCol(1:2)/Box%size(1:2)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(Box%wave, this%wave_norm, xCol(3), exp_kzCol_tab)

        mColOverL(:) = 2._DP*PI * mCol(1:2)/Box%size(1:2)

        do ky = 0, Box%wave(2)

            wave_vector(2) = real(ky, DP)

            do kx = -Box_wave1_sym(Box%wave, ky, 0), Box%wave(1)

                wave_vector(1) = real(kx, DP)

                kMcol_z = this%wave_norm(kx, ky) * mCol(3)
                k_dot_mCol = dot_product(wave_vector, mColOverL)
                
                exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky)
                exp_kzCol = exp_kzCol_tab(abs(kx), ky)
                
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                                              cmplx(+kMcol_z, k_dot_mCol, DP) * &
                                              exp_IkxCol * cmplx(exp_kzCol, 0._DP, DP)
                
                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                                               cmplx(-kMcol_z, k_dot_mCol, DP) * &
                                               exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)

            end do
            
        end do

    end subroutine Electronic_Layer_Correction_update_structures_exchange

end module class_dipolarSpheres
