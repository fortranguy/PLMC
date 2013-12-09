!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, real_zero, consist_tiny
use data_constants, only : PI
use data_box, only : Ndim, Lsize, Kmax, out_permittivity
use data_particles, only : dipol_Ncol
use data_monteCarlo, only : dipol_move_delta, dipol_move_rejectFix, dipol_rotate_delta, &
                            dipol_rotate_deltaMax, dipol_rotate_rejectFix, dipol_Nwidom, &
                            dipol_reInit_iStep
use data_potential, only : dipol_rCut, dipol_dr, dipol_alpha
use data_neighbourCells, only : NnearCell
use data_distribution, only : dipol_snap_factor
use module_physics, only : set_discrete_length, distVec_PBC, Kmax1_sym, Kmax2_sym, fourier_i
use class_neighbourCells
use class_hardSpheres

implicit none

private

    type, extends(HardSpheres), public :: DipolarSpheres

        private
        
        ! Particles
        
        real(DP), dimension(:, :), allocatable, public :: orientations !< dipolar orientations 
                                                                       !< of all particles
        
        ! Monte-Carlo
        real(DP), public :: rotate_delta !< rotation
        real(DP) :: rotate_deltaSave
        real(DP) :: rotate_deltaMax
        real(DP) :: rotate_rejectFix
        integer :: reInit_iStep

        ! Potential
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :), allocatable :: Epot_real_tab !< tabulation : real short-range
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: Epot_reci_weight
        integer :: NwaveVectors
        complex(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: structure
        real(DP), dimension(Ndim) :: totalMoment
        
    contains

        !> Construction and destruction of the class
        procedure :: init_particles => DipolarSpheres_init_particles
        procedure :: init_changes => DipolarSpheres_init_changes
        procedure :: construct => DipolarSpheres_construct
        procedure :: destroy => DipolarSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: print_report => DipolarSpheres_print_report
        
        !> Take a snap shot of the configuration : orientations
        procedure :: snap_orientations => DipolarSpheres_snap_orientations
        
        !> Adapt the rotation rotate_delta during thermalisation
        procedure :: adapt_rotate_delta => DipolarSpheres_adapt_rotate_delta
        procedure :: set_rotate_delta => DipolarSpheres_set_rotate_delta
        procedure :: get_rotate_delta => DipolarSpheres_get_rotate_delta
        
        procedure :: get_reInit_iStep => DipolarSpheres_get_reInit_iStep
        
        !> Potential energy
        !>     Real
        procedure, private :: Epot_real_true => DipolarSpheres_Epot_real_true
        procedure, private :: Epot_real_set_parameters => DipolarSpheres_Epot_real_set_parameters
        procedure, private :: Epot_real_set_tab => DipolarSpheres_Epot_real_set_tab
        procedure, private :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure, private :: Epot_real_interpol => DipolarSpheres_Epot_real_interpol
        procedure, private :: Epot_real_pair => DipolarSpheres_Epot_real_pair
        procedure :: Epot_real_solo => DipolarSpheres_Epot_real_solo
        procedure, private :: Epot_real => DipolarSpheres_Epot_real
        !>     Reciprocal : init
        procedure, private :: Epot_reci_init_weight => DipolarSpheres_Epot_reci_init_weight
        procedure, private :: Epot_reci_init_structure => DipolarSpheres_Epot_reci_init_structure
        procedure, private :: Epot_reci_init => DipolarSpheres_Epot_reci_init
        procedure, private :: Epot_reci_get_structure_modulus => &
                              DipolarSpheres_Epot_reci_get_structure_modulus
        procedure :: Epot_reci_reInit_structure => DipolarSpheres_Epot_reci_reInit_structure
        procedure :: Epot_reci_count_waveVectors => DipolarSpheres_Epot_reci_count_waveVectors
        !>     Reciprocal : delta
        procedure :: deltaEpot_reci_move => DipolarSpheres_deltaEpot_reci_move
        procedure :: reci_update_structure_move => &
                              DipolarSpheres_reci_update_structure_move
        procedure :: deltaEpot_reci_rotate => DipolarSpheres_deltaEpot_reci_rotate
        procedure :: reci_update_structure_rotate => &
                              DipolarSpheres_reci_update_structure_rotate
        procedure :: deltaEpot_reci_exchange => DipolarSpheres_deltaEpot_reci_exchange
        !>     Reciprocal : total
        procedure, private :: Epot_reci => DipolarSpheres_Epot_reci
        !>     Self
        procedure :: Epot_self_solo => DipolarSpheres_Epot_self_solo
        procedure, private :: Epot_self => DipolarSpheres_Epot_self
        !>     Total moment
        procedure :: init_totalMoment => DipolarSpheres_init_totalMoment
        procedure :: reInit_totalMoment => DipolarSpheres_reInit_totalMoment
        procedure :: update_totalMoment_rotate => DipolarSpheres_update_totalMoment_rotate
        !>     Boundary conditions
        procedure :: deltaEpot_bound_rotate => DipolarSpheres_deltaEpot_bound_rotate
        procedure :: deltaEpot_bound_exchange => DipolarSpheres_deltaEpot_bound_exchange
        procedure, private :: Epot_bound => DipolarSpheres_Epot_bound
        !>     Total
        procedure :: init_potential => DipolarSpheres_init_potential
        procedure :: Epot_conf => DipolarSpheres_Epot_conf
        procedure :: test_consist => DipolarSpheres_test_consist
        
    end type DipolarSpheres
    
contains

    subroutine DipolarSpheres_init_particles(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        this%rMin = 1._DP
        this%radius = this%rMin/2._DP
        this%Ncol = dipol_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        allocate(this%orientations(Ndim, this%Ncol))
    
    end subroutine DipolarSpheres_init_particles
    
    subroutine DipolarSpheres_init_changes(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        this%move_delta = dipol_move_delta
        this%move_deltaSave = this%move_delta
        this%move_rejectFix = dipol_move_rejectFix
        
        this%rotate_delta = dipol_rotate_delta
        this%rotate_deltaSave = this%rotate_delta
        this%rotate_deltaMax = dipol_rotate_deltaMax
        this%rotate_rejectFix = dipol_rotate_rejectFix
        
        this%Nwidom = dipol_Nwidom
        
    end subroutine DipolarSpheres_init_changes

    subroutine DipolarSpheres_construct(this)
    
        class(DipolarSpheres), intent(out) :: this
        
        real(DP), dimension(Ndim) :: cell_size
        
        this%name = "dipol"        
        write(output_unit, *) this%name, " class construction"
    
        call this%init_particles()
        
        ! Snapshot
        this%snap_factor = dipol_snap_factor
        
        call this%init_changes()
        
        ! Neighbour Cells
        cell_size(:) = this%rMin
        call this%sameCells%construct(cell_size, this%rCut) !< same kind
    
    end subroutine DipolarSpheres_construct
    
    subroutine DipolarSpheres_destroy(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        call this%HardSpheres%destroy()
        
        if (allocated(this%orientations)) then
            deallocate(this%orientations)
        end if
        
        if (allocated(this%Epot_real_tab)) then
            deallocate(this%Epot_real_tab)
        end if
    
    end subroutine DipolarSpheres_destroy
    
    !> Report
    
    subroutine DipolarSpheres_print_report(this, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        call this%HardSpheres%print_report(report_unit)

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    dr = ", this%dr
        write(report_unit, *) "    reInit_iStep = ", this%reInit_iStep
        write(report_unit, *) "    NwaveVectors = ", this%NwaveVectors
        
    end subroutine DipolarSpheres_print_report
    
    !> Configuration state : orientations
      
    subroutine DipolarSpheres_snap_orientations(this, iStep, snap_unit)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%orientations(:, iCol)
            end do
            
        end if

    end subroutine DipolarSpheres_snap_orientations
    
    !> Adaptation of rotate_delta during the thermalisation
    
    pure subroutine DipolarSpheres_adapt_rotate_delta(this, reject)
    
        class(DipolarSpheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: rotate_delta_eps = 0.05_DP
        real(DP), parameter :: rotate_reject_eps = 0.1_DP * rotate_delta_eps
        real(DP), parameter :: more = 1._DP+rotate_delta_eps
        real(DP), parameter :: less = 1._DP-rotate_delta_eps

        
        if (reject < this%rotate_rejectFix - rotate_reject_eps) then
        
            this%rotate_delta = this%rotate_delta * more
  
            if (this%rotate_delta > this%rotate_deltaMax) then
                this%rotate_delta = this%rotate_deltaMax
            end if
            
        else if (reject > this%rotate_rejectFix + rotate_reject_eps) then
        
            this%rotate_delta = this%rotate_delta * less
            
        end if
    
    end subroutine DipolarSpheres_adapt_rotate_delta
    
    subroutine DipolarSpheres_set_rotate_delta(this, reject, report_unit)
    
        class(DipolarSpheres), intent(inout) :: this    
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        if (reject < real_zero) then
            write(error_unit, *) this%name, " :    Warning : rotate_delta adaptation problem."
            this%rotate_delta = this%rotate_deltaSave
            write(error_unit, *) "default rotate_delta :", this%rotate_delta
        end if
        
        if (this%rotate_delta > this%rotate_deltaMax) then
            write(error_unit, *) this%name, " :   Warning : rotate_delta too big."
            this%rotate_delta = this%rotate_deltaMax
            write(error_unit, *) "big rotate_delta :", this%rotate_delta
        end if
        
        write(report_unit, *) "Rotation :"
        write(report_unit, *) "    rotate_delta = ", this%rotate_delta
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%rotate_rejectFix)/this%rotate_rejectFix
    
    end subroutine DipolarSpheres_set_rotate_delta
    
    !> Accessor : rotate_delta
    
    pure function DipolarSpheres_get_rotate_delta(this) result(get_rotate_delta)
        
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: get_rotate_delta
        
        get_rotate_delta = this%rotate_delta
        
    end function DipolarSpheres_get_rotate_delta
    
    !> Accessor : reInit_iStep
    
    pure function DipolarSpheres_get_reInit_iStep(this) result (get_reInit_iStep)
    
        class(DipolarSpheres), intent(in) :: this
        integer :: get_reInit_iStep
    
        get_reInit_iStep = this%reInit_iStep
        
    end function DipolarSpheres_get_reInit_iStep

    ! Real : short-range interaction ---------------------------------------------------------------

    !> \f[ B(r) = \frac{\mathrm{erfc}(\alpha r)}{r^3} +
    !>           2\frac{\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 r^2}}{r^2} \f]
    !> \f[ C(r) = 3\frac{\mathrm{erfc}(\alpha r)}{r^5} +
    !>            2\frac{\alpha}{\sqrt{\pi}}(2\alpha^2 + \frac{3}{r^2})
    !>                                     \frac{e^{-\alpha^2 r^2}}{r^2} \f]

    pure function DipolarSpheres_Epot_real_true(this, r) result(Epot_real_true)

        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP), dimension(2) :: Epot_real_true

        real(DP) :: alpha

        alpha = this%alpha

        Epot_real_true(1) = erfc(alpha*r)/r**3 + 2._DP*alpha/sqrt(PI) * exp(-alpha**2*r**2) / r**2

        Epot_real_true(2) = 3._DP*erfc(alpha*r)/r**5 + &
                            2._DP*alpha/sqrt(PI) * (2._DP*alpha**2+3._DP/r**2) * &
                            exp(-alpha**2*r**2) / r**2

    end function DipolarSpheres_Epot_real_true
    
    !> Initialisation : look-up (tabulation) table
    
    pure subroutine DipolarSpheres_Epot_real_set_tab(this)
    
        class(DipolarSpheres), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
        real(DP) :: alpha
        
        alpha = this%alpha
       
        ! cut
        do i = this%iMin, this%iCut
        
            r_i = real(i, DP)*this%dr
            
            this%Epot_real_tab(i, :) = this%Epot_real_true(r_i)
            
        end do
        
        ! shift        
        this%Epot_real_tab(:, 1) = this%Epot_real_tab(:, 1) - this%Epot_real_tab(this%iCut, 1)
        this%Epot_real_tab(:, 2) = this%Epot_real_tab(:, 2) - this%Epot_real_tab(this%iCut, 2)

    end subroutine DipolarSpheres_Epot_real_set_tab
    
    !> Initialisation
    
    subroutine DipolarSpheres_Epot_real_set_parameters(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        this%rCut = dipol_rCut
        this%dr = dipol_dr
        call set_discrete_length(this%rMin, this%dr)
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr) + 1
        
    end subroutine DipolarSpheres_Epot_real_set_parameters
    
    subroutine DipolarSpheres_Epot_real_init(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        call this%Epot_real_set_parameters()
        
        allocate(this%Epot_real_tab(this%iMin:this%iCut, 2))
        call this%Epot_real_set_tab()
    
    end subroutine DipolarSpheres_Epot_real_init

    !> Print the tabulated values
    
    subroutine DipolarSpheres_Epot_real_print(this, Epot_unit)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        integer :: i
        real(DP) :: r_i

        do i = this%iMin, this%iCut
            r_i = real(i, DP)*this%dr
            write(Epot_unit, *) r_i, this%Epot_real_tab(i, :)
        end do

    end subroutine DipolarSpheres_Epot_real_print
    
    !> Linear interpolation

    pure function DipolarSpheres_Epot_real_interpol(this, r) result(Epot_real_interpol)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: r
        real(DP), dimension(2) :: Epot_real_interpol
        
        integer :: i
        real(DP) :: r_i
       
        if (r < this%rCut) then
       
            i = int(r/this%dr)
            r_i = real(i, DP)*this%dr
            Epot_real_interpol(:) = this%Epot_real_tab(i, :) + (r-r_i)/this%dr * &
                                   (this%Epot_real_tab(i+1, :) - this%Epot_real_tab(i, :))
           
        else
       
            Epot_real_interpol(:) = 0._DP
           
        end if
        
    end function DipolarSpheres_Epot_real_interpol

    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B(r_{ij}) - 
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C(r_{ij}) \f]
    
    pure function DipolarSpheres_Epot_real_pair(this, mCol_i, mCol_j, rVec_ij, r_ij) &
                  result(Epot_real_pair)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol_i, mCol_j
        real(DP), dimension(:), intent(in) :: rVec_ij
        real(DP), intent(in) :: r_ij
        real(DP) :: Epot_real_pair
        
        real(DP), dimension(2) :: Epot_coeff
        
        Epot_coeff(1) = dot_product(mCol_i, mCol_j)
        Epot_coeff(2) =-dot_product(mCol_i, rVec_ij) * dot_product(mCol_j, rVec_ij)
        
        Epot_real_pair = dot_product(Epot_coeff, this%Epot_real_interpol(r_ij))
    
    end function DipolarSpheres_Epot_real_pair
    
    !> Energy of 1 dipole with others
    
    pure function DipolarSpheres_Epot_real_solo(this, iCol, xCol, mCol) result(Epot_real_solo)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol
        real(DP), dimension(:), intent(in) :: xCol, mCol
        real(DP) :: Epot_real_solo

        integer :: jCol
        real(DP), dimension(Ndim) :: rVec_ij
        real(DP) :: r_ij

        Epot_real_solo = 0._DP

        do jCol = 1, this%Ncol

            if (jCol /= iCol) then

                rVec_ij = distVec_PBC(xCol(:), this%positions(:, jCol))
                r_ij = norm2(rVec_ij)

                Epot_real_solo = Epot_real_solo + &
                                 this%Epot_real_pair(mCol, this%orientations(:, jCol), rVec_ij, r_ij)

            end if

        end do
        
    end function DipolarSpheres_Epot_real_solo
    
    !> Total real energy
    
    pure function DipolarSpheres_Epot_real(this) result(Epot_real)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_real
        
        integer :: iCol, jCol
        real(DP), dimension(Ndim) :: rVec_ij
        real(DP) :: r_ij
    
        Epot_real = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    rVec_ij = distVec_PBC(this%positions(:, iCol), this%positions(:, jCol))
                    r_ij = norm2(rVec_ij)
                    
                    Epot_real = Epot_real + &
                                this%Epot_real_pair(this%orientations(:, iCol), &
                                                    this%orientations(:, jCol), rVec_ij, r_ij)
                    
                end if
            end do
        end do
        
        Epot_real = 0.5_DP*Epot_real
    
    end function DipolarSpheres_Epot_real

    ! Reciprocal : long-range interaction ----------------------------------------------------------
    
    !> \f[ 
    !>      w(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}{\alpha^2} \sum_d \frac{k_d^2}{L_d}}}
    !>                                {\sum_d \frac{k_d^2}{L_d}}
    !> \f]
    
    pure subroutine DipolarSpheres_Epot_reci_init_weight(this)
        
        class(DipolarSpheres), intent(inout) :: this
        
        integer :: kx, ky, kz
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: kOverL

        do kz = -Kmax(3), Kmax(3)
            waveVector(3) = real(kz, DP)
        
        do ky = -Kmax(2), Kmax(2)
            waveVector(2) = real(ky, DP)
        
        do kx = -Kmax(1), Kmax(1)
            waveVector(1) = real(kx, DP)

            if (kx**2 + ky**2 + kz**2 /= 0) then
            
                kOverL = norm2(waveVector(:)/Lsize(:))

                this%Epot_reci_weight(kx, ky, kz) = exp(-PI**2/this%alpha**2 * kOverL**2) / kOverL**2

            else

                this%Epot_reci_weight(kx, ky, kz) = 0._DP

            end if

        end do
            
        end do
        
        end do
        
    end subroutine DipolarSpheres_Epot_reci_init_weight
    
    !> Structure factor init :
    !> \f[ 
    !>      S(\vec{k}) = \sum_{i} (\vec{k}\cdot\vec{\mu}_i) e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f]
    !> We will also use a restricted definition later :
    !> \f[ 
    !>      S_l(\vec{k}) = \sum_{i \neq l} (\vec{k}\cdot\vec{\mu}_i) e^{+i\vec{k}\cdot\vec{x}_i}
    !> \f].

    pure subroutine DipolarSpheres_Epot_reci_init_structure(this)

        class(DipolarSpheres), intent(inout) :: this

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3

        real(DP), dimension(Ndim) :: xColOverL, mColOverL
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        integer :: iCol

        this%structure(:, :, :) = cmplx(0._DP, 0._DP, DP)

        do iCol = 1, this%Ncol
        
            xColOverL(:) = 2._DP*PI * this%positions(:, iCol)/Lsize(:)            
            call fourier_i(Kmax(1), xColOverL(1), exp_Ikx_1)
            call fourier_i(Kmax(2), xColOverL(2), exp_Ikx_2)
            call fourier_i(Kmax(3), xColOverL(3), exp_Ikx_3)
            
            mColOverL(:) = this%orientations(:, iCol)/Lsize(:)
        
            do kz = -Kmax(3), Kmax(3)

                waveVector(3) = real(kz, DP)

            do ky = -Kmax(2), Kmax(2)

                waveVector(2) = real(ky, DP)

            do kx = -Kmax(1), Kmax(1)

                waveVector(1) = real(kx, DP)
                
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                
                k_dot_mCol = dot_product(waveVector, mColOverL)
                          
                this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                             cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine DipolarSpheres_Epot_reci_init_structure
    
    subroutine DipolarSpheres_Epot_reci_init(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        call this%Epot_reci_init_weight()
        call this%Epot_reci_init_structure()
    
    end subroutine DipolarSpheres_Epot_reci_init
    
    !> To calculate the drift of the strucutre factor

    pure function DipolarSpheres_Epot_reci_get_structure_modulus(this) &
                  result(Epot_reci_get_structure_modulus)

        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci_get_structure_modulus

        integer :: kx, ky, kz

        Epot_reci_get_structure_modulus = 0._DP

        do kz = 0, Kmax(3)
            do ky = -Kmax2_sym(kz), Kmax(2)
                do kx = -Kmax1_sym(ky, kz), Kmax(1)
                
                    Epot_reci_get_structure_modulus = Epot_reci_get_structure_modulus + &
                                                      abs(this%structure(kx, ky, kz))

                end do
            end do
        end do

    end function DipolarSpheres_Epot_reci_get_structure_modulus
    
    !> Reinitialise the structure factor and print the drift
    
    subroutine DipolarSpheres_Epot_reci_reInit_structure(this, iStep, modulus_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reInit
        
        modulus_drifted = this%Epot_reci_get_structure_modulus()
        call this%Epot_reci_init_structure()
        modulus_reInit = this%Epot_reci_get_structure_modulus()
        
        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)
    
    end subroutine DipolarSpheres_Epot_reci_reInit_structure

    ! Count the number of wave vectors

    subroutine DipolarSpheres_Epot_reci_count_waveVectors(this, waveVectors_unit)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: waveVectors_unit
        
        integer :: kx, ky, kz

        this%NwaveVectors = 0

        do kz = 0, Kmax(3)

            do ky = -Kmax2_sym(kz), Kmax(2)

                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    if (kx**2 + ky**2 + kz**2 /= 0) then

                        write(waveVectors_unit, *) kx, ky, kz
                        write(waveVectors_unit, *)
                        write(waveVectors_unit, *)

                        this%NwaveVectors = this%NwaveVectors + 1

                    end if

                end do

            end do

        end do

    end subroutine DipolarSpheres_Epot_reci_count_waveVectors

    !> Move

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta M^2
    !> w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta M^2 = 2\Re[
    !>                  (\vec{\mu}_l\cdot\vec{k})
    !>                  (e^{-i\vec{k}\cdot\vec{x}^\prime_l} - e^{-i\vec{k}\cdot\vec{x}_l})
    !>                  S_l(\vec{k})
    !>               ]
    !> \f]

    pure function DipolarSpheres_deltaEpot_reci_move(this, xOld, xNew, mCol) &
                  result(deltaEpot_reci_move)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_reci_move
        
        real(DP), dimension(Ndim) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz

        xNewOverL(:) = 2._DP*PI * xNew(:)/Lsize(:)
        call fourier_i(Kmax(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Kmax(2), xNewOverL(2), exp_IkxNew_2)
        call fourier_i(Kmax(3), xNewOverL(3), exp_IkxNew_3)
        
        xOldOverL(:) = 2._DP*PI * xOld(:)/Lsize(:)
        call fourier_i(Kmax(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Kmax(2), xOldOverL(2), exp_IkxOld_2)
        call fourier_i(Kmax(3), xOldOverL(3), exp_IkxOld_3)

        mColOverL(:) = mCol(:)/Lsize(:)

        deltaEpot_reci_move = 0._DP

        do kz = 0, Kmax(3) ! symmetry : half wave vectors -> double Energy

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2) 

                waveVector(2) = real(ky, DP)
            
                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)
                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)

                    realPart = k_dot_mCol * real((conjg(exp_IkxNew) - conjg(exp_IkxOld)) * &
                    (this%structure(kx, ky, kz) - cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxOld), DP)

                    deltaEpot_reci_move = deltaEpot_reci_move + &
                                          2._DP * this%Epot_reci_weight(kx, ky, kz) * realPart

                end do
            
            end do
        
        end do

        deltaEpot_reci_move = 4._DP*PI/product(Lsize) * deltaEpot_reci_move

    end function DipolarSpheres_deltaEpot_reci_move

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    pure subroutine DipolarSpheres_reci_update_structure_move(this, xOld, xNew, mCol)

        class(DipolarSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        
        real(DP), dimension(Ndim) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz

        xNewOverL(:) = 2._DP*PI * xNew(:)/Lsize(:)
        call fourier_i(Kmax(1), xNewOverL(1), exp_IkxNew_1)
        call fourier_i(Kmax(2), xNewOverL(2), exp_IkxNew_2)
        call fourier_i(Kmax(3), xNewOverL(3), exp_IkxNew_3)
        
        xOldOverL(:) = 2._DP*PI * xOld(:)/Lsize(:)
        call fourier_i(Kmax(1), xOldOverL(1), exp_IkxOld_1)
        call fourier_i(Kmax(2), xOldOverL(2), exp_IkxOld_2)
        call fourier_i(Kmax(3), xOldOverL(3), exp_IkxOld_3)

        mColOverL(:) = mCol(:)/Lsize(:)

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)
                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                                                          
                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                        cmplx(k_dot_mCol, 0._DP, DP) * (exp_IkxNew - exp_IkxOld)

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_reci_update_structure_move
    
    !> Rotate

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta M^2
    !>                                       w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta M^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2\Re\{
    !>                  [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>                  e^{-i \vec{k} \cdot \vec{x}_l} S_l(\vec{k})
    !>               \}
    !> \f]
    
    pure function DipolarSpheres_deltaEpot_reci_rotate(this, xCol, mOld, mNew) &
                  result(deltaEpot_reci_rotate)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_reci_rotate

        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol

        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mNew, k_dot_mOld
        integer :: kx, ky, kz

        xColOverL(:) = 2._DP*PI * xCol(:)/Lsize(:)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Kmax(3), xColOverL(3), exp_IkxCol_3)

        mNewOverL(:) = mNew(:)/Lsize(:)
        mOldOverL(:) = mOld(:)/Lsize(:)

        deltaEpot_reci_rotate = 0._DP

        do kz = 0, Kmax(3) ! symmetry : half wave vectors -> double Energy

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)
            
                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)

                    k_dot_mNew = dot_product(waveVector, mNewOverL)
                    k_dot_mOld = dot_product(waveVector, mOldOverL)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    realPart = k_dot_mNew**2 - k_dot_mOld**2
                    realPart = realPart + 2._DP * (k_dot_mNew - k_dot_mOld) * real(conjg(exp_IkxCol) * &
                    (this%structure(kx, ky, kz) - k_dot_mOld * exp_IkxCol), DP)

                    deltaEpot_reci_rotate = deltaEpot_reci_rotate + &
                                            this%Epot_reci_weight(kx, ky, kz) * realPart

                end do

            end do

        end do

        deltaEpot_reci_rotate = 4._DP*PI/product(Lsize) * deltaEpot_reci_rotate

    end function DipolarSpheres_deltaEpot_reci_rotate

    !> Update moment -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = [\vec{k}\cdot(\vec{\mu}_l^\prime - \vec{\mu}_l)] 
    !>                          e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    pure subroutine DipolarSpheres_reci_update_structure_rotate(this, xCol, mOld, mNew)

        class(DipolarSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew

        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol

        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_deltaMcol
        integer :: kx, ky, kz

        xColOverL(:) = 2._DP*PI * xCol(:)/Lsize(:)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Kmax(3), xColOverL(3), exp_IkxCol_3)

        mNewOverL(:) = mNew(:)/Lsize(:)
        mOldOverL(:) = mOld(:)/Lsize(:)

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)
                    k_dot_deltaMcol = dot_product(waveVector, mNewOverL - mOldOverL)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                        cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_IkxCol

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_reci_update_structure_rotate

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
    
    !> Summary : only the sign of \f[\vec{\mu}\f] changes.

    pure function DipolarSpheres_deltaEpot_reci_exchange(this, xCol, mCol) &
                  result(deltaEpot_reci_exchange)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_reci_exchange
        
        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mColOverL
        
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol
        
        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        
        xColOverL(:) = 2._DP*PI * xCol(:)/Lsize(:)
        call fourier_i(Kmax(1), xColOverL(1), exp_IkxCol_1)
        call fourier_i(Kmax(2), xColOverL(2), exp_IkxCol_2)
        call fourier_i(Kmax(3), xColOverL(3), exp_IkxCol_3)
        
        mColOverL(:) = mCol(:)/Lsize(:)
        
        deltaEpot_reci_exchange = 0._DP
        
        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)
            
                do kx = -Kmax1_sym(ky, kz), Kmax(1)
                
                    waveVector(1) = real(kx, DP)                    
                    k_dot_mCol = dot_product(waveVector, mColOverL)
                                                  
                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)
                    
                    realPart = k_dot_mCol * (k_dot_mCol + 2._DP * &
                    real(this%structure(kx, ky, kz) * conjg(exp_IkxCol), DP))

                    deltaEpot_reci_exchange = deltaEpot_reci_exchange + &
                                              this%Epot_reci_weight(kx, ky, kz) * realPart
                   
                end do
            
            end do
        
        end do
        
        deltaEpot_reci_exchange = 4._DP*PI/product(Lsize) * deltaEpot_reci_exchange

    end function DipolarSpheres_deltaEpot_reci_exchange
    
    !> Total reciprocal energy
    
    pure function DipolarSpheres_Epot_reci(this) result(Epot_reci)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci

        integer :: kx, ky, kz

        Epot_reci = 0._DP

        do kz = -Kmax(3), Kmax(3)

            do ky = -Kmax(2), Kmax(2)

                do kx = -Kmax(1), Kmax(1)

                    Epot_reci = Epot_reci + this%Epot_reci_weight(kx, ky, kz) * &
                                real(this%structure(kx, ky, kz)*conjg(this%structure(kx, ky, kz)), DP)
        
                end do

            end do

        end do
        
        Epot_reci = 2._DP*PI/product(Lsize) * Epot_reci
        
    end function DipolarSpheres_Epot_reci
    
    ! Self : correction ----------------------------------------------------------------------------
    
    !> Self energy of 1 dipole
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function DipolarSpheres_Epot_self_solo(this, mCol) result(Epot_self_solo)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: Epot_self_solo
        
        Epot_self_solo = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(mCol, mCol)
    
    end function DipolarSpheres_Epot_self_solo
    
    !> Total self energy
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    pure function DipolarSpheres_Epot_self(this) result(Epot_self)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_self

        integer :: iCol        
        
        Epot_self = 0._DP
        do iCol = 1, this%Ncol
            Epot_self = Epot_self + this%Epot_self_solo(this%orientations(:, iCol))
        end do
        
    end function DipolarSpheres_Epot_self

    ! Total moment ---------------------------------------------------------------------------------

    !> Total dipole moment :
    !> \f[ \vec{M} = \sum_j \vec{\mu}_j \f]
    !> \f[ \vec{M}_l = \sum_{j \neq l} \vec{\mu}_j \f]
    
    pure subroutine DipolarSpheres_init_totalMoment(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        integer :: iCol
        
        this%totalMoment(:) = 0._DP
        
        do iCol = 1, this%Ncol
        
            this%totalMoment(:) = this%totalMoment(:) + this%orientations(:, iCol)
        
        end do        
    
    end subroutine DipolarSpheres_init_totalMoment

    !> Reinitialise the total moment factor and print the drift

    subroutine DipolarSpheres_reInit_totalMoment(this, iStep, modulus_unit)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit

        real(DP) :: modulus_drifted, modulus_reInit

        modulus_drifted = norm2(this%totalMoment(:))
        call this%init_totalMoment()
        modulus_reInit = norm2(this%totalMoment(:))

        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)

    end subroutine DipolarSpheres_reInit_totalMoment

    !> Rotation

    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]

    pure subroutine DipolarSpheres_update_totalMoment_rotate(this, mOld, mNew)

        class(DipolarSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mOld, mNew

        this%totalMoment(:) = this%totalMoment(:) + mNew(:) - mOld(:)

    end subroutine DipolarSpheres_update_totalMoment_rotate

    ! Boundary conditions : shape-dependent -------------------------------------------------------
    
    !> Exchange
    
    !> Difference of Energy : add
    !> \f[
    !>      \Delta U_{N \rightarrow N+1} = \frac{2\pi}{(2\epsilon_s+1)V} [
    !>                                         (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          2\vec{\mu}) \cdot \vec{M}^N
    !>                                     ]
    !> \f]
    
    !> Difference of Energy : remove
    !> \f[
    !>      \Delta U_{N \rightarrow N-1} = \frac{2\pi}{(2\epsilon_s+1)V} [
    !>                                          (\vec{\mu}_{N+1} \cdot \vec{\mu}_{N+1})
    !>                                          -2\vec{\mu}) \cdot \vec{M}^N
    !>                                      ]
    !> \f]
    
    pure function DipolarSpheres_deltaEpot_bound_exchange(this, mCol) &
                  result (deltaEpot_bound_exchange)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_bound_exchange
        
        deltaEpot_bound_exchange = dot_product(mCol, mCol) + &
                                   2._DP * dot_product(mCol, this%totalMoment)
                          
        deltaEpot_bound_exchange = 2._DP*PI / (2._DP*out_permittivity+1._DP) / product(Lsize) * &
                                   deltaEpot_bound_exchange
    
    end function DipolarSpheres_deltaEpot_bound_exchange
    
    !> Rotation
    
    !> Difference of Energy
    !> \f[
    !>      \Delta U = \frac{2\pi}{(2\epsilon_s+1)V} [
    !>                      (\vec{\mu}^\prime_l \cdot \vec{\mu}^\prime_l) -
    !>                      (\vec{\mu}_l \cdot \vec{\mu}_l) +
    !>                      2 (\vec{\mu}^\prime_l - \vec{\mu}_l) \cdot \vec{M}_l
    !>                 ]
    !> \f]
    
    pure function DipolarSpheres_deltaEpot_bound_rotate(this, mOld, mNew) &
                  result (deltaEpot_bound_rotate)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_bound_rotate
        
        deltaEpot_bound_rotate = dot_product(mNew, mNew) - dot_product(mOld, mOld) + &
                                 2._DP*dot_product(mNew-mOld, this%totalMoment-mOld)
                          
        deltaEpot_bound_rotate = 2._DP*PI / (2._DP*out_permittivity+1._DP) / product(Lsize) * &
                                 deltaEpot_bound_rotate
    
    end function DipolarSpheres_deltaEpot_bound_rotate
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{(2\epsilon_s+1)V} | \vec{M}|^2
    !> \f]
    
    pure function DipolarSpheres_Epot_bound(this) result(Epot_bound)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_bound
        
        Epot_bound = 2._DP * PI / (2._DP*out_permittivity + 1._DP) / product(Lsize) * &
                     dot_product(this%totalMoment, this%totalMoment)
    
    end function DipolarSpheres_Epot_bound
    
    !> Potential energy initialisation
    
    subroutine DipolarSpheres_init_potential(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        this%alpha = dipol_alpha
        
        call this%Epot_real_init()
        
        call this%Epot_reci_init()
        this%reInit_iStep = dipol_reInit_iStep
        
    end subroutine DipolarSpheres_init_potential

    !> Total potential energy
    
    pure function DipolarSpheres_Epot_conf(this) result(Epot_conf)
    
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
        
        Epot_conf = this%Epot_real() + this%Epot_reci() - this%Epot_self() + this%Epot_bound()
    
    end function DipolarSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine DipolarSpheres_test_consist(this, Epot, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
        real(DP) :: difference
        
        Epot_conf = this%Epot_conf()
        difference = abs((Epot_conf-Epot)/Epot_conf)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    relative difference = ", difference
        
        if (difference > consist_tiny) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine DipolarSpheres_test_consist

end module class_dipolarSpheres
