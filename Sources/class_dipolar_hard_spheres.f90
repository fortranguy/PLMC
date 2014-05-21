!> \brief Description of the Dipolar_Hard_Spheres class

module class_dipolar_hard_spheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero, consist_tiny
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use data_monte_carlo, only: dipol_rotate_delta, dipol_rotate_deltaMax, dipol_rotate_rejectFix
use data_potential, only: dipol_rMin_factor, dipol_real_rCut_factor, dipol_real_dr, dipol_alpha_factor
use data_neighbour_cells, only: NnearCell
use data_distribution, only: snap_ratio
use module_types_micro, only: Box_Dimensions, Particle_Index
use module_physics_micro, only: set_discrete_length, PBC_vector, Box_wave1_sym, Box_wave2_sym, &
                                fourier_i, exchange_sign
use module_data, only: test_data_found
use class_hard_spheres

implicit none

private

    type, extends(Hard_Spheres), public :: Dipolar_Hard_Spheres

        private
        
        ! Particles
        real(DP), dimension(:, :), allocatable, public :: all_orientations
        ! Potential
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :, :), allocatable :: Epot_reci_weight
        integer :: NwaveVectors
        complex(DP), dimension(:, :, :), allocatable :: Epot_reci_structure
        real(DP), dimension(Ndim) :: totalMoment
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Dipolar_Hard_Spheres_construct
        procedure, private :: set_particles => Dipolar_Hard_Spheres_set_particles
        procedure :: destroy => Dipolar_Hard_Spheres_destroy
        
        !> Write a report of the component in a file
        procedure :: write_report => Dipolar_Hard_Spheres_write_report
        
        !> Accessor & Mutator
        procedure :: get_orientation => Dipolar_Hard_Spheres_get_orientation
        procedure :: set_orientation => Dipolar_Hard_Spheres_set_orientation
        procedure :: Epot_set_alpha => Dipolar_Hard_Spheres_Epot_set_alpha
        
        !> Take a snap shot of the configuration: orientations
        procedure :: write_snap_orientations => Dipolar_Hard_Spheres_write_snap_orientations
        
        !> Potential energy
        !>     Real
        !>     Reciprocal: init
        procedure, private :: set_Epot_reci_weight => Dipolar_Hard_Spheres_set_Epot_reci_weight
        procedure, private :: set_Epot_reci_structure => Dipolar_Hard_Spheres_set_Epot_reci_structure
        procedure, private :: set_Epot_reci => Dipolar_Hard_Spheres_set_Epot_reci
        procedure, private :: Epot_reci_get_structure_modulus => &
                              Dipolar_Hard_Spheres_Epot_reci_get_structure_modulus
        procedure :: reset_Epot_reci_structure => Dipolar_Hard_Spheres_reset_Epot_reci_structure
        procedure :: Epot_reci_count_waveVectors => Dipolar_Hard_Spheres_Epot_reci_count_waveVectors
        !>     Reciprocal: delta
        procedure :: deltaEpot_reci_move => Dipolar_Hard_Spheres_deltaEpot_reci_move
        procedure :: reci_update_structure_move => &
                              Dipolar_Hard_Spheres_reci_update_structure_move
        procedure :: deltaEpot_reci_rotate => Dipolar_Hard_Spheres_deltaEpot_reci_rotate
        procedure :: reci_update_structure_rotate => &
                              Dipolar_Hard_Spheres_reci_update_structure_rotate
        procedure :: deltaEpot_reci_exchange => Dipolar_Hard_Spheres_deltaEpot_reci_exchange
        !>     Reciprocal: total
        procedure, private :: Epot_reci => Dipolar_Hard_Spheres_Epot_reci
        !>     Self
        procedure :: Epot_self_solo => Dipolar_Hard_Spheres_Epot_self_solo
        procedure, private :: Epot_self => Dipolar_Hard_Spheres_Epot_self
        !>     Total moment
        procedure, private :: set_totalMoment => Dipolar_Hard_Spheres_set_totalMoment
        procedure :: reset_totalMoment => Dipolar_Hard_Spheres_reset_totalMoment
        procedure :: update_totalMoment_rotate => Dipolar_Hard_Spheres_update_totalMoment_rotate
        !>     Boundary conditions
        procedure :: deltaEpot_bound_rotate => Dipolar_Hard_Spheres_deltaEpot_bound_rotate
        procedure :: deltaEpot_bound_exchange => Dipolar_Hard_Spheres_deltaEpot_bound_exchange
        procedure, private :: Epot_bound => Dipolar_Hard_Spheres_Epot_bound
        !>     Total
        procedure :: set_Epot => Dipolar_Hard_Spheres_set_Epot
        procedure :: Epot_conf => Dipolar_Hard_Spheres_Epot_conf
        
    end type Dipolar_Hard_Spheres
    
contains

    subroutine Dipolar_Hard_Spheres_construct(this, json)
    
        class(Dipolar_Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        this%name = "dipol"
        write(output_unit, *) this%name, " class construction"
    
        call this%set_particles(json)        
        call this%Hard_Spheres%set_snap(json)
    
    end subroutine Dipolar_Hard_Spheres_construct

    subroutine Dipolar_Hard_Spheres_set_particles(this, json)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        this%diameter = 1._DP ! = u_length
        
        data_name = "Particles.Dipoles.number of particles"
        call json%get(data_name, this%num_particles, found)
        call test_data_found(data_name, found)
        allocate(this%all_positions(Ndim, this%num_particles))
        allocate(this%all_orientations(Ndim, this%num_particles))
        
        data_name = "Particles.Dipoles.number of Widom particles"
        call json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
        
    end subroutine Dipolar_Hard_Spheres_set_particles
    
    subroutine Dipolar_Hard_Spheres_destroy(this)    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        
        call this%Hard_Spheres%destroy()
        if (allocated(this%all_orientations)) deallocate(this%all_orientations)
        if (allocated(this%Epot_reci_weight)) deallocate(this%Epot_reci_weight)
        if (allocated(this%Epot_reci_structure)) deallocate(this%Epot_reci_structure)    
    end subroutine Dipolar_Hard_Spheres_destroy
    
    !> Report
    
    subroutine Dipolar_Hard_Spheres_write_report(this, report_unit)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        call this%Hard_Spheres%write_report(report_unit)

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    NwaveVectors = ", this%NwaveVectors
        
    end subroutine Dipolar_Hard_Spheres_write_report
    
    !> Accessors & Mutators
    
    pure function Dipolar_Hard_Spheres_get_orientation(this, i_particle) result(get_orientation)
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(Ndim) :: get_orientation
        
        get_orientation(:) = this%all_orientations(:, i_particle)
    end function Dipolar_Hard_Spheres_get_orientation
    
    subroutine Dipolar_Hard_Spheres_set_orientation(this, i_particle, orientation)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(:), intent(in) :: orientation
        
        this%all_orientations(:, i_particle) = orientation(:)
    end subroutine Dipolar_Hard_Spheres_set_orientation
    
    subroutine Dipolar_Hard_Spheres_Epot_set_alpha(this, Box_size)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        
        this%alpha = dipol_alpha_factor / Box_size(1)
    end subroutine Dipolar_Hard_Spheres_Epot_set_alpha
    
    !> Configuration state: orientations
      
    subroutine Dipolar_Hard_Spheres_write_snap_orientations(this, iStep, snap_unit)
        
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(iStep, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_orientations(:, i_particle)
            end do
        end if

    end subroutine Dipolar_Hard_Spheres_write_snap_orientations    
    
    ! Reciprocal: long-range interaction ----------------------------------------------------------
    
    !> \f[
    !>      w(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}{\alpha^2} \sum_{d=1}^3 \frac{k_d^2}{L_d}}}
    !>                                {\sum_{d=1}^3 \frac{k_d^2}{L_d}}
    !> \f]
    
    pure subroutine Dipolar_Hard_Spheres_set_Epot_reci_weight(this, Box)
        
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        
        integer :: kx, ky, kz
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: kOverL

        do kz = -Box%wave(3), Box%wave(3)
            waveVector(3) = real(kz, DP)
        
        do ky = -Box%wave(2), Box%wave(2)
            waveVector(2) = real(ky, DP)
        
        do kx = -Box%wave(1), Box%wave(1)
            waveVector(1) = real(kx, DP)

            if (kx**2 + ky**2 + kz**2 /= 0) then
                kOverL = norm2(waveVector(:)/Box%size(:))
                this%Epot_reci_weight(kx, ky, kz) = exp(-PI**2/this%alpha**2 * kOverL**2) / kOverL**2
            else
                this%Epot_reci_weight(kx, ky, kz) = 0._DP
            end if

        end do
            
        end do
        
        end do
        
    end subroutine Dipolar_Hard_Spheres_set_Epot_reci_weight
    
    !> Structure factor init :
    !> \f[
    !>      S(\vec{k}) = \sum_{i} (\vec{k}\cdot\vec{\mu}_i) e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f]
    !> We will also use a restricted definition later :
    !> \f[
    !>      S_\underline{l}(\vec{k}) = \sum_{i \neq l} (\vec{k}\cdot\vec{\mu}_i)
    !>                                 e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f].

    pure subroutine Dipolar_Hard_Spheres_set_Epot_reci_structure(this, Box)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3

        real(DP), dimension(Ndim) :: xColOverL, mColOverL
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        integer :: i_particle

        this%Epot_reci_structure(:, :, :) = cmplx(0._DP, 0._DP, DP)

        do i_particle = 1, this%num_particles
        
            xColOverL(:) = 2._DP*PI * this%all_positions(:, i_particle)/Box%size(:)
            call fourier_i(Box%wave(1), xColOverL(1), exp_Ikx_1)
            call fourier_i(Box%wave(2), xColOverL(2), exp_Ikx_2)
            call fourier_i(Box%wave(3), xColOverL(3), exp_Ikx_3)
            
            mColOverL(:) = this%all_orientations(:, i_particle)/Box%size(:)
        
            do kz = -Box%wave(3), Box%wave(3)
                waveVector(3) = real(kz, DP)

            do ky = -Box%wave(2), Box%wave(2)
                waveVector(2) = real(ky, DP)
                
            do kx = -Box%wave(1), Box%wave(1)
                waveVector(1) = real(kx, DP)
                              
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                k_dot_mCol = dot_product(waveVector, mColOverL)
                          
                this%Epot_reci_structure(kx, ky, kz) = this%Epot_reci_structure(kx, ky, kz) + &
                                                       cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine Dipolar_Hard_Spheres_set_Epot_reci_structure
    
    pure subroutine Dipolar_Hard_Spheres_set_Epot_reci(this, Box)
        
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        
        if (allocated(this%Epot_reci_weight)) deallocate(this%Epot_reci_weight)
        allocate(this%Epot_reci_weight(-Box%wave(1):Box%wave(1), -Box%wave(2):Box%wave(2), &
                                       -Box%wave(3):Box%wave(3)))
        
        if (allocated(this%Epot_reci_structure)) deallocate(this%Epot_reci_structure)
        allocate(this%Epot_reci_structure(-Box%wave(1):Box%wave(1), -Box%wave(2):Box%wave(2), &
                                          -Box%wave(3):Box%wave(3)))
        
        call this%set_Epot_reci_weight(Box)
        call this%set_Epot_reci_structure(Box)
    
    end subroutine Dipolar_Hard_Spheres_set_Epot_reci
    
    !> To calculate the drift of the strucutre factor

    pure function Dipolar_Hard_Spheres_Epot_reci_get_structure_modulus(this, Box_wave) &
                  result(Epot_reci_get_structure_modulus)

        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, dimension(:), intent(in) :: Box_wave
        real(DP) :: Epot_reci_get_structure_modulus

        integer :: kx, ky, kz

        Epot_reci_get_structure_modulus = 0._DP

        do kz = 0, Box_wave(3)
            do ky = -Box_wave2_sym(Box_wave, kz), Box_wave(2)
                do kx = -Box_wave1_sym(Box_wave, ky, kz), Box_wave(1)
                    Epot_reci_get_structure_modulus = Epot_reci_get_structure_modulus + &
                                                      abs(this%Epot_reci_structure(kx, ky, kz))
                end do
            end do
        end do

    end function Dipolar_Hard_Spheres_Epot_reci_get_structure_modulus
    
    !> Reinitialise the structure factor and write the drift
    
    subroutine Dipolar_Hard_Spheres_reset_Epot_reci_structure(this, Box, iStep, modulus_unit)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reInit
        
        modulus_drifted = this%Epot_reci_get_structure_modulus(Box%wave)
        call this%set_Epot_reci_structure(Box)
        modulus_reInit = this%Epot_reci_get_structure_modulus(Box%wave)
        
        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)
    
    end subroutine Dipolar_Hard_Spheres_reset_Epot_reci_structure

    ! Count the number of wave vectors

    subroutine Dipolar_Hard_Spheres_Epot_reci_count_waveVectors(this, Box_wave, waveVectors_unit)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, dimension(:), intent(in) :: Box_wave
        integer, intent(in) :: waveVectors_unit
        
        integer :: kx, ky, kz

        this%NwaveVectors = 0

        do kz = 0, Box_wave(3)
            do ky = -Box_wave2_sym(Box_wave, kz), Box_wave(2)
                do kx = -Box_wave1_sym(Box_wave, ky, kz), Box_wave(1)
                    if (kx**2 + ky**2 + kz**2 /= 0) then

                        write(waveVectors_unit, *) kx, ky, kz
                        write(waveVectors_unit, *)
                        write(waveVectors_unit, *)

                        this%NwaveVectors = this%NwaveVectors + 1

                    end if
                end do
            end do
        end do

    end subroutine Dipolar_Hard_Spheres_Epot_reci_count_waveVectors

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

    pure function Dipolar_Hard_Spheres_deltaEpot_reci_move(this, Box, old, new) &
                  result(deltaEpot_reci_move)

        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: deltaEpot_reci_move
        
        real(DP), dimension(Ndim) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim) :: mColOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz

        xNewOverL(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), xNewOverL(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), xNewOverL(3), exp_IkxNew_3)
        
        xOldOverL(:) = 2._DP*PI * old%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), xOldOverL(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), xOldOverL(3), exp_IkxOld_3)

        mColOverL(:) = new%orientation(:)/Box%size(:)

        deltaEpot_reci_move = 0._DP

        do kz = 0, Box%wave(3) ! symmetry: half wave vectors -> double Energy
            waveVector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                waveVector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)

                    realPart = k_dot_mCol * real((conjg(exp_IkxNew) - conjg(exp_IkxOld)) * &
                        (this%Epot_reci_structure(kx, ky, kz) - cmplx(k_dot_mCol, 0._DP, DP) * &
                        exp_IkxOld), DP)

                    deltaEpot_reci_move = deltaEpot_reci_move + &
                                          2._DP * this%Epot_reci_weight(kx, ky, kz) * realPart

                end do
            
            end do
        
        end do

        deltaEpot_reci_move = 4._DP*PI/product(Box%size) * deltaEpot_reci_move

    end function Dipolar_Hard_Spheres_deltaEpot_reci_move

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    pure subroutine Dipolar_Hard_Spheres_reci_update_structure_move(this, Box, old, new)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        
        real(DP), dimension(Ndim) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim) :: mColOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz

        xNewOverL(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Box%wave(2), xNewOverL(2), exp_IkxNew_2)
        call fourier_i(Box%wave(3), xNewOverL(3), exp_IkxNew_3)
        
        xOldOverL(:) = 2._DP*PI * old%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Box%wave(2), xOldOverL(2), exp_IkxOld_2)
        call fourier_i(Box%wave(3), xOldOverL(3), exp_IkxOld_3)

        mColOverL(:) = new%orientation(:)/Box%size(:)

        do kz = 0, Box%wave(3)
            waveVector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                waveVector(2) = real(ky, DP)

                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_mCol = dot_product(waveVector, mColOverL)
                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                                                          
                    this%Epot_reci_structure(kx, ky, kz) = this%Epot_reci_structure(kx, ky, kz) + &
                        cmplx(k_dot_mCol, 0._DP, DP) * (exp_IkxNew - exp_IkxOld)

                end do
                
            end do
            
        end do

    end subroutine Dipolar_Hard_Spheres_reci_update_structure_move
    
    !> Rotate

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta S^2
    !>                                       w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta S^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2\Re\{
    !>                  [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>                  e^{-i \vec{k} \cdot \vec{x}_l} S_\underline{l}(\vec{k})
    !>               \}
    !> \f]
    
    pure function Dipolar_Hard_Spheres_deltaEpot_reci_rotate(this, Box, old, new) &
                  result(deltaEpot_reci_rotate)

        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: deltaEpot_reci_rotate

        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol

        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mNew, k_dot_mOld
        integer :: kx, ky, kz

        xColOverL(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Box%wave(3), xColOverL(3), exp_IkxCol_3)

        mNewOverL(:) = new%orientation(:)/Box%size(:)
        mOldOverL(:) = old%orientation(:)/Box%size(:)

        deltaEpot_reci_rotate = 0._DP

        do kz = 0, Box%wave(3) ! symmetry: half wave vectors -> double Energy
            waveVector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                waveVector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    waveVector(1) = real(kx, DP)

                    k_dot_mNew = dot_product(waveVector, mNewOverL)
                    k_dot_mOld = dot_product(waveVector, mOldOverL)
                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    realPart = k_dot_mNew**2 - k_dot_mOld**2
                    realPart = realPart + 2._DP * (k_dot_mNew-k_dot_mOld) * real(conjg(exp_IkxCol) * &
                        (this%Epot_reci_structure(kx, ky, kz) - k_dot_mOld * exp_IkxCol), DP)

                    deltaEpot_reci_rotate = deltaEpot_reci_rotate + &
                                            this%Epot_reci_weight(kx, ky, kz) * realPart

                end do

            end do

        end do

        deltaEpot_reci_rotate = 4._DP*PI/product(Box%size) * deltaEpot_reci_rotate

    end function Dipolar_Hard_Spheres_deltaEpot_reci_rotate

    !> Update moment -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = [\vec{k}\cdot(\vec{\mu}_l^\prime - \vec{\mu}_l)]
    !>                          e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    pure subroutine Dipolar_Hard_Spheres_reci_update_structure_rotate(this, Box, old, new)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new

        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol

        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_deltaMcol
        integer :: kx, ky, kz

        xColOverL(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Box%wave(3), xColOverL(3), exp_IkxCol_3)

        mNewOverL(:) = new%orientation(:)/Box%size(:)
        mOldOverL(:) = old%orientation(:)/Box%size(:)

        do kz = 0, Box%wave(3)
            waveVector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                waveVector(2) = real(ky, DP)

                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_deltaMcol = dot_product(waveVector, mNewOverL - mOldOverL)
                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    this%Epot_reci_structure(kx, ky, kz) = this%Epot_reci_structure(kx, ky, kz) + &
                        cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_IkxCol

                end do
                
            end do
            
        end do

    end subroutine Dipolar_Hard_Spheres_reci_update_structure_rotate

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

    pure function Dipolar_Hard_Spheres_deltaEpot_reci_exchange(this, Box, particle) &
                  result(deltaEpot_reci_exchange)

        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP) :: deltaEpot_reci_exchange
        
        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mColOverL
        
        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol
        
        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        
        xColOverL(:) = 2._DP*PI * particle%position(:)/Box%size(:)
        call fourier_i(Box%wave(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Box%wave(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Box%wave(3), xColOverL(3), exp_IkxCol_3)
        
        mColOverL(:) = exchange_sign(particle%add) * particle%orientation(:)/Box%size(:)
        
        deltaEpot_reci_exchange = 0._DP
        
        do kz = 0, Box%wave(3)
            waveVector(3) = real(kz, DP)

            do ky = -Box_wave2_sym(Box%wave, kz), Box%wave(2)
                waveVector(2) = real(ky, DP)
            
                do kx = -Box_wave1_sym(Box%wave, ky, kz), Box%wave(1)
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_mCol = dot_product(waveVector, mColOverL)
                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)
                    
                    realPart = k_dot_mCol * (k_dot_mCol + 2._DP * &
                    real(this%Epot_reci_structure(kx, ky, kz) * conjg(exp_IkxCol), DP))

                    deltaEpot_reci_exchange = deltaEpot_reci_exchange + &
                                              this%Epot_reci_weight(kx, ky, kz) * realPart
                   
                end do
            
            end do
        
        end do
        
        deltaEpot_reci_exchange = 4._DP*PI/product(Box%size) * deltaEpot_reci_exchange

    end function Dipolar_Hard_Spheres_deltaEpot_reci_exchange
    
    !> Total reciprocal energy
    
    pure function Dipolar_Hard_Spheres_Epot_reci(this, Box) result(Epot_reci)
        
        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP) :: Epot_reci

        integer :: kx, ky, kz

        Epot_reci = 0._DP

        do kz = -Box%wave(3), Box%wave(3)
            do ky = -Box%wave(2), Box%wave(2)
                do kx = -Box%wave(1), Box%wave(1)
                    Epot_reci = Epot_reci + this%Epot_reci_weight(kx, ky, kz) * &
                                            real(this%Epot_reci_structure(kx, ky, kz) * &
                                            conjg(this%Epot_reci_structure(kx, ky, kz)), DP)
                end do
            end do
        end do
        
        Epot_reci = 2._DP*PI/product(Box%size) * Epot_reci
        
    end function Dipolar_Hard_Spheres_Epot_reci
    
    ! Self: correction ----------------------------------------------------------------------------
    
    !> Self energy of 1 dipole
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Dipolar_Hard_Spheres_Epot_self_solo(this, mCol) result(Epot_self_solo)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: Epot_self_solo
        
        Epot_self_solo = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(mCol, mCol)
    
    end function Dipolar_Hard_Spheres_Epot_self_solo
    
    !> Total self energy
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function Dipolar_Hard_Spheres_Epot_self(this) result(Epot_self)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP) :: Epot_self

        integer :: i_particle
        
        Epot_self = 0._DP
        do i_particle = 1, this%num_particles
            Epot_self = Epot_self + this%Epot_self_solo(this%all_orientations(:, i_particle))
        end do
        
    end function Dipolar_Hard_Spheres_Epot_self

    ! Total moment ---------------------------------------------------------------------------------

    !> Total dipole moment :
    !> \f[ \vec{M} = \sum_j \vec{\mu}_j \f]
    !> \f[ \vec{M}_\underline{l} = \sum_{j \neq l} \vec{\mu}_j \f]
    
    pure subroutine Dipolar_Hard_Spheres_set_totalMoment(this)
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer :: i_particle
        this%totalMoment(:) = 0._DP
        do i_particle = 1, this%num_particles
            this%totalMoment(:) = this%totalMoment(:) + this%all_orientations(:, i_particle)
        end do
    end subroutine Dipolar_Hard_Spheres_set_totalMoment

    !> Reinitialise the total moment factor and write the drift

    subroutine Dipolar_Hard_Spheres_reset_totalMoment(this, iStep, modulus_unit)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reInit

        modulus_drifted = norm2(this%totalMoment(:))
        call this%set_totalMoment()
        modulus_reInit = norm2(this%totalMoment(:))

        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)

    end subroutine Dipolar_Hard_Spheres_reset_totalMoment

    !> Rotation

    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]

    pure subroutine Dipolar_Hard_Spheres_update_totalMoment_rotate(this, mOld, mNew)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mOld, mNew

        this%totalMoment(:) = this%totalMoment(:) + mNew(:) - mOld(:)

    end subroutine Dipolar_Hard_Spheres_update_totalMoment_rotate

    ! Boundary conditions: shape-dependent -------------------------------------------------------
    
    !> Exchange
    
    !> Difference of Energy: add
    !> \f[
    !>      \Delta U_{N \rightarrow N+1} = \frac{2\pi}{3V} [
    !>                                         (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          +2(\vec{\mu} \cdot \vec{M}_N)
    !>                                     ]
    !> \f]
    
    !> Difference of Energy: remove
    !> \f[
    !>      \Delta U_{N \rightarrow N-1} = \frac{2\pi}{3V} [
    !>                                          (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          -2(\vec{\mu} \cdot \vec{M}_N)
    !>                                      ]
    !> \f]
    
    pure function Dipolar_Hard_Spheres_deltaEpot_bound_exchange(this, Box_size, mCol) &
                  result (deltaEpot_bound_exchange)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_bound_exchange
        
        deltaEpot_bound_exchange = dot_product(mCol, mCol) + &
                                   2._DP * dot_product(mCol, this%totalMoment)
                          
        deltaEpot_bound_exchange = 2._DP*PI/3._DP/product(Box_size) * deltaEpot_bound_exchange
    
    end function Dipolar_Hard_Spheres_deltaEpot_bound_exchange
    
    !> Rotation
    
    !> Difference of Energy
    !> \f[
    !>      \Delta U = \frac{2\pi}{3V} [
    !>                      (\vec{\mu}^\prime_l \cdot \vec{\mu}^\prime_l) -
    !>                      (\vec{\mu}_l \cdot \vec{\mu}_l) +
    !>                      2 (\vec{\mu}^\prime_l - \vec{\mu}_l) \cdot \vec{M}_\underline{l}
    !>                 ]
    !> \f]
    
    pure function Dipolar_Hard_Spheres_deltaEpot_bound_rotate(this, Box_size, mOld, mNew) &
                  result (deltaEpot_bound_rotate)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_bound_rotate
        
        deltaEpot_bound_rotate = dot_product(mNew, mNew) - dot_product(mOld, mOld) + &
                                 2._DP*dot_product(mNew-mOld, this%totalMoment-mOld)
                          
        deltaEpot_bound_rotate = 2._DP*PI/3._DP/product(Box_size) * deltaEpot_bound_rotate
    
    end function Dipolar_Hard_Spheres_deltaEpot_bound_rotate
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{3V} | \vec{M}|^2
    !> \f]
    
    pure function Dipolar_Hard_Spheres_Epot_bound(this, Box_size) result(Epot_bound)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP) :: Epot_bound
        
        Epot_bound = 2._DP*PI/3._DP/product(Box_size) * dot_product(this%totalMoment, this%totalMoment)
    
    end function Dipolar_Hard_Spheres_Epot_bound
    
    ! Total potential energy ----------------------------------------------------------------------
    
    !> Potential energy initialisation
    
    subroutine Dipolar_Hard_Spheres_set_Epot(this, Box)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        
        call this%Epot_set_alpha(Box%size)
        call this%set_Epot_reci(Box)
        call this%set_totalMoment()
        
    end subroutine Dipolar_Hard_Spheres_set_Epot

    !> Total potential energy of a configuration
    
    pure function Dipolar_Hard_Spheres_Epot_conf(this, Box) result(Epot_conf)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP) :: Epot_conf
        
        Epot_conf = this%Epot_reci(Box) - this%Epot_self() + &
                    this%Epot_bound(Box%size) ! this%Epot_real(Box%size)
    
    end function Dipolar_Hard_Spheres_Epot_conf

end module class_dipolar_hard_spheres
