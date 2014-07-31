!> \brief Description of the Electronic Layer Correction class

module class_dipolarSpheres

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use data_box, only: num_dimensions
use module_types_micro, only: Box_Parameters, Particle_Index
use module_physics_micro, only: fourier_i, Box_wave1_sym, Box_wave2_sym, exchange_sign, set_exp_kz
use class_hard_spheres, only: Dipolar_Hard_Spheres

implicit none

private

    type, public :: Electronic_Layer_Correction

        private
        
        integer :: NwaveVectors
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2)) :: norm_k
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2)) :: weight
        complex(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2)) :: structure_plus, structure_minus
        
    contains

        !>     ELC: init
        procedure, private :: init_norm_k => Electronic_Layer_Correction_init_norm_k
        procedure, private :: init_weight => Electronic_Layer_Correction_init_weight
        procedure, private :: init_structures => Electronic_Layer_Correction_init_structures
        procedure, private :: init => Electronic_Layer_Correction_init
        procedure, private :: get_structures_modulus => &
                              Electronic_Layer_Correction_get_structures_modulus
        procedure :: reInit_structures => Electronic_Layer_Correction_reInit_structures
        procedure :: count_waveVectors => Electronic_Layer_Correction_count_waveVectors
        !>     ELC: delta
        procedure :: deltamove => Electronic_Layer_Correction_deltamove
        procedure :: update_structures_move => &
                     Electronic_Layer_Correction_update_structures_move
        procedure :: deltarotate => Electronic_Layer_Correction_deltarotate
        procedure :: update_structures_rotate => &
                     Electronic_Layer_Correction_update_structures_rotate
        procedure :: deltaexchange => Electronic_Layer_Correction_deltaexchange
        !>     ELC: total
        procedure, private :: ELC => Electronic_Layer_Correction_ELC
        
    end type Electronic_Layer_Correction
    
contains
    
    !> Initialisation of the ``structure factors''
    
    pure subroutine Electronic_Layer_Correction_init_norm_k(this)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
    
        integer :: kx, ky
        real(DP), dimension(Ndim-1) :: waveVector
        
        do ky = -Kmax(2), Kmax(2)
            waveVector(2) = real(ky, DP)
        
        do kx = -Kmax(1), Kmax(1)
            waveVector(1) = real(kx, DP)
            
            this%norm_k(kx, ky) = 2._DP*PI * norm2(waveVector(:)/Lsize(1:Ndim-1))
            
        end do
        
        end do
    
    end subroutine Electronic_Layer_Correction_init_norm_k
    
    !> \f[
    !>      w(\vec{k}^{2D}) = \frac{1}{k^{2D}(e^{k^{2D}L_z} - 1)}
    !> \f]
    
    pure subroutine Electronic_Layer_Correction_init_weight(this)
        
        class(Electronic_Layer_Correction), intent(inout) :: this
        
        integer :: kx, ky
        
        do ky = -Kmax(2), Kmax(2)
        
        do kx = -Kmax(1), Kmax(1)

            if (kx**2 + ky**2 /= 0) then

                this%weight(kx, ky) = 1._DP / this%norm_k(kx, ky) / &
                                               (exp(this%norm_k(kx, ky)*Lsize(3)) - 1._DP)

            else

                this%weight(kx, ky) = 0._DP

            end if

        end do
            
        end do
        
    end subroutine Electronic_Layer_Correction_init_weight
    
    !> Structure factor init :
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

    pure subroutine Electronic_Layer_Correction_init_structures(this)

        class(Electronic_Layer_Correction), intent(inout) :: this

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        real(DP) :: exp_kzCol

        real(DP), dimension(Ndim-1) :: xColOverL
        real(DP), dimension(Ndim-1) :: mColOverL
        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky
        integer :: iCol

        this%structure_plus(:, :) = cmplx(0._DP, 0._DP, DP)
        this%structure_minus(:, :) = cmplx(0._DP, 0._DP, DP)

        do iCol = 1, this%Ncol
        
            xColOverL(:) = 2._DP*PI * this%positions(1:2, iCol)/Lsize(1:2)
            call fourier_i(Kmax(1), xColOverL(1), exp_Ikx_1)
            call fourier_i(Kmax(2), xColOverL(2), exp_Ikx_2)
            
            mColOverL(:) = 2._DP*PI * this%orientations(1:2, iCol)/Lsize(1:2)
        
            do ky = -Kmax(2), Kmax(2)

                waveVector(2) = real(ky, DP)

            do kx = -Kmax(1), Kmax(1)

                waveVector(1) = real(kx, DP)
            
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky)
                exp_kzCol = exp(this%norm_k(kx, ky) * this%positions(3, iCol))

                k_dot_mCol = dot_product(waveVector, mColOverL)
                kMcol_z = this%norm_k(kx, ky) * this%orientations(3, iCol)
                          
                this%structure_plus(kx, ky) = this%structure_plus(kx, ky) + &
                                              cmplx(+kMcol_z, k_dot_mCol, DP) * &
                                              exp_IkxCol * cmplx(exp_kzCol, 0._DP, DP)

                this%structure_minus(kx, ky) = this%structure_minus(kx, ky) + &
                                               cmplx(-kMcol_z, k_dot_mCol, DP) * &
                                               exp_IkxCol / cmplx(exp_kzCol, 0._DP, DP)
            
            end do
            
            end do
            
        end do

    end subroutine Electronic_Layer_Correction_init_structures
    
    !> Initialisation
    
    subroutine Electronic_Layer_Correction_init(this)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
    
        call this%init_norm_k()
        call this%init_weight()
        call this%init_structures()
    
    end subroutine Electronic_Layer_Correction_init
    
    !> To calculate the drift of the structure factor

    pure function Electronic_Layer_Correction_get_structures_modulus(this) &
                  result(get_structures_modulus)

        class(Electronic_Layer_Correction), intent(in) :: this
        real(DP), dimension(2) :: get_structures_modulus

        integer :: kx, ky

        get_structures_modulus(:) = 0._DP

        do ky = 0, Kmax(2)
            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                get_structures_modulus(1) = get_structures_modulus(1) + &
                                                     abs(this%structure_plus(kx, ky))
                get_structures_modulus(2) = get_structures_modulus(2) + &
                                                     abs(this%structure_minus(kx, ky))
            end do
        end do

    end function Electronic_Layer_Correction_get_structures_modulus
    
    !> Reinitialise the structure factor and print the drift
    
    subroutine Electronic_Layer_Correction_reInit_structures(this, iStep, modulus_unit)
    
        class(Electronic_Layer_Correction), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP), dimension(2) :: modulus_drifted, modulus_reInit
        
        modulus_drifted(:) = this%get_structures_modulus()
        call this%init_structures()
        modulus_reInit(:) = this%get_structures_modulus()
        
        write(modulus_unit, *) iStep, abs(modulus_reInit(:) - modulus_drifted(:))
    
    end subroutine Electronic_Layer_Correction_reInit_structures
    
    ! Count the number of wave vectors

    subroutine Electronic_Layer_Correction_count_waveVectors(this, waveVectors_unit)

        class(Electronic_Layer_Correction), intent(inout) :: this
        integer, intent(in) :: waveVectors_unit
        
        integer :: kx, ky

        this%NwaveVectors = 0

        do ky = 0, Kmax(2)
            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                if (kx**2 + ky**2 /= 0) then
                    write(waveVectors_unit, *) kx, ky
                    write(waveVectors_unit, *)
                    write(waveVectors_unit, *)

                    this%NwaveVectors = this%NwaveVectors + 1

                end if
            end do
        end do

    end subroutine Electronic_Layer_Correction_count_waveVectors
    
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

    pure function Electronic_Layer_Correction_deltamove(this, xOld, xNew, mCol) &
                  result(deltamove)

        class(Electronic_Layer_Correction), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltamove

        real(DP), dimension(Ndim-1) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim-1) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP) :: exp_IkxNew
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzNew_tab
        real(DP) :: exp_kzNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP) :: exp_IkxOld
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzOld_tab
        real(DP) :: exp_kzOld

        complex(DP) :: structure_i, delta_structure_i
        
        real(DP) :: realPart1, realPart2
        
        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky

        xNewOverL(:) = 2._DP*PI * xNew(1:Ndim-1)/Lsize(1:Ndim-1)
        call fourier_i(Kmax(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Kmax(2), xNewOverL(2), exp_IkxNew_2)
        call set_exp_kz(this%norm_k, xNew(3), exp_kzNew_tab)
        
        xOldOverL(:) = 2._DP*PI * xOld(1:Ndim-1)/Lsize(1:Ndim-1)
        call fourier_i(Kmax(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Kmax(2), xOldOverL(2), exp_IkxOld_2)
        call set_exp_kz(this%norm_k, xOld(3), exp_kzOld_tab)

        mColOverL(:) = 2._DP*PI * mCol(1:Ndim-1)/Lsize(1:Ndim-1)

        deltamove = 0._DP

        do ky = 0, Kmax(2)
            waveVector(2) = real(ky, DP)
        
            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                waveVector(1) = real(kx, DP)
                
                kMcol_z = this%norm_k(kx, ky) * mCol(3)
                k_dot_mCol = dot_product(waveVector, mColOverL)
                
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
                
                deltamove = deltamove + this%weight(kx, ky) * 2._DP * &
                                     (realPart1 + realPart2)

            end do
        
        end do

        deltamove = 2._DP * PI/product(Lsize(1:2)) * deltamove

    end function Electronic_Layer_Correction_deltamove

    !> Update position -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = (\pm k^{2D}\mu_{l ,z} + i(\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))
    !>                       (e^{\pm k^{2D}z_l^\prime} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D\prime}_l)} -
    !>                        e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)})
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structures_move(this, xOld, xNew, mCol)

        class(Electronic_Layer_Correction), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        
        real(DP), dimension(Ndim-1) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim-1) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP) :: exp_IkxNew
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzNew_tab
        real(DP) :: exp_kzNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP) :: exp_IkxOld
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzOld_tab
        real(DP) :: exp_kzOld

        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: k_dot_mCol, kMcol_z
        integer :: kx, ky

        xNewOverL(:) = 2._DP*PI * xNew(1:Ndim-1)/Lsize(1:Ndim-1)
        call fourier_i(Kmax(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Kmax(2), xNewOverL(2), exp_IkxNew_2)
        call set_exp_kz(this%norm_k, xNew(3), exp_kzNew_tab)
        
        xOldOverL(:) = 2._DP*PI * xOld(1:Ndim-1)/Lsize(1:Ndim-1)
        call fourier_i(Kmax(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Kmax(2), xOldOverL(2), exp_IkxOld_2)
        call set_exp_kz(this%norm_k, xOld(3), exp_kzOld_tab)

        mColOverL(:) = 2._DP*PI * mCol(1:Ndim-1)/Lsize(1:Ndim-1)

        do ky = 0, Kmax(2)
            waveVector(2) = real(ky, DP)

            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                waveVector(1) = real(kx, DP)

                k_dot_mCol = dot_product(waveVector, mColOverL)
                kMcol_z = this%norm_k(kx, ky) * mCol(3)

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
    
    pure function Electronic_Layer_Correction_deltarotate(this, xCol, mOld, mNew) &
                  result(deltarotate)

        class(Electronic_Layer_Correction), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltarotate

        real(DP), dimension(Ndim-1) :: xColOverL
        real(DP), dimension(Ndim-1) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol

        complex(DP) :: structure_i, delta_structure_i
        real(DP) :: realPart0, realPart1, realPart2
        
        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: kMnew_z, kMold_z
        real(DP) :: k_dot_mNew, k_dot_mOld
        integer :: kx, ky

        xColOverL(:) = 2._DP*PI * xCol(1:2)/Lsize(1:2)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(this%norm_k, xCol(3), exp_kzCol_tab)

        mNewOverL(:) = 2._DP*PI * mNew(1:2)/Lsize(1:2)
        mOldOverL(:) = 2._DP*PI * mOld(1:2)/Lsize(1:2)

        deltarotate = 0._DP

        do ky = 0, Kmax(2)
            waveVector(2) = real(ky, DP)
        
            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                waveVector(1) = real(kx, DP)

                kMnew_z = this%norm_k(kx, ky) * mNew(3)
                k_dot_mNew = dot_product(waveVector, mNewOverL)
                
                kMold_z = this%norm_k(kx, ky) * mOld(3)
                k_dot_mOld = dot_product(waveVector, mOldOverL)

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

                deltarotate = deltarotate + this%weight(kx, ky) * 2._DP * &
                                     (realPart0 + realPart1 + realPart2)

            end do

        end do

        deltarotate = 2._DP * PI/product(Lsize(1:2)) * deltarotate

    end function Electronic_Layer_Correction_deltarotate

    !> Update moment -> update the ``structure factors''
    !>  \f[
    !>      \Delta S_{\pm} = [\pm (k^{2D}\mu_{l ,z}^\prime - k^{2D}\mu_{l ,z}) +
    !>                       i((\vec{k}^{2D}\cdot\vec{\mu}^{2D\prime}_l) -
    !>                         (\vec{k}^{2D}\cdot\vec{\mu}^{2D}_l))]
    !>                       e^{\pm k^{2D}z_l} e^{i(\vec{k}^{2D}\cdot\vec{x}^{2D}_l)}
    !>  \f]

    pure subroutine Electronic_Layer_Correction_update_structures_rotate(this, xCol, mOld, mNew)

        class(Electronic_Layer_Correction), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew

        real(DP), dimension(Ndim-1) :: xColOverL
        real(DP), dimension(Ndim-1) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol

        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: kdeltaMcol_z, k_dot_deltaMcol
        integer :: kx, ky

        xColOverL(:) = 2._DP*PI * xCol(1:Ndim-1)/Lsize(1:Ndim-1)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(this%norm_k, xCol(3), exp_kzCol_tab)

        mNewOverL(:) = 2._DP*PI * mNew(1:Ndim-1)/Lsize(1:Ndim-1)
        mOldOverL(:) = 2._DP*PI * mOld(1:Ndim-1)/Lsize(1:Ndim-1)

        do ky = 0, Kmax(2)
            waveVector(2) = real(ky, DP)

            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                waveVector(1) = real(kx, DP)

                exp_kzCol = exp_kzCol_tab(abs(kx), ky)
                exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky)

                kdeltaMcol_z = this%norm_k(kx, ky) * (mNew(3) - mOld(3))
                k_dot_deltaMcol = dot_product(waveVector, mNewOverL - mOldOverL)
                
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

    pure function Electronic_Layer_Correction_deltaexchange(this, xCol, mCol) &
                  result(deltaexchange)

        class(Electronic_Layer_Correction), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaexchange
        
        real(DP), dimension(Ndim-1) :: xColOverL
        real(DP), dimension(Ndim-1) :: mColOverL
        
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP) :: exp_IkxCol
        real(DP), dimension(0:Kmax(1), 0:Kmax(2)) :: exp_kzCol_tab
        real(DP) :: exp_kzCol
        
        complex(DP) :: structure_i
        real(DP) :: realPart0, realPart1, realPart2
        
        real(DP), dimension(Ndim-1) :: waveVector
        real(DP) :: kMcol_z, k_dot_mCol
        integer :: kx, ky
        
        xColOverL(:) = 2._DP*PI * xCol(1:2)/Lsize(1:2)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call set_exp_kz(this%norm_k, xCol(3), exp_kzCol_tab)
        
        mColOverL(:) = 2._DP*PI * mCol(1:2)/Lsize(1:2)
        
        deltaexchange = 0._DP

        do ky = 0, Kmax(2)
            waveVector(2) = real(ky, DP)
        
            do kx = -Kmax1_sym(ky, 0), Kmax(1)
                waveVector(1) = real(kx, DP)
                
                kMcol_z = this%norm_k(kx, ky) * mCol(3)
                k_dot_mCol = dot_product(waveVector, mColOverL)
                                                
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

                deltaexchange = deltaexchange + this%weight(kx, ky) * &
                                         2._DP * (realPart0 + realPart1 + realPart2)
                
            end do
        
        end do

        deltaexchange = 2._DP * PI/product(Lsize(1:2)) * deltaexchange

    end function Electronic_Layer_Correction_deltaexchange
    
    !> Total ELC energy
    
    pure function Electronic_Layer_Correction_ELC(this) result(ELC)
        
        class(Electronic_Layer_Correction), intent(in) :: this
        real(DP) :: ELC

        complex(DP) :: structure_product
        integer kx, ky

        ELC = 0._DP

        do ky = -Kmax(2), Kmax(2)
            do kx = -Kmax(1), Kmax(1)
                structure_product = this%structure_plus(kx, ky)*conjg(this%structure_minus(kx, ky))
                ELC = ELC + this%weight(kx, ky) * 2._DP*real(structure_product, DP)
            end do
        end do
        
        ELC = PI/product(Lsize(1:2)) * ELC
        
    end function Electronic_Layer_Correction_ELC

end module class_dipolarSpheres
