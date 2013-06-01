!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding, only : C_int, C_double
use data_constants
use data_cell
use data_particles
use data_potentiel
use data_mc
use data_neighbours
use data_distrib
use mod_physics
use class_neighbours
use class_mixingPotential
use class_spheres
use mod_ewald_reci

implicit none

private

    type, extends(Spheres), public :: DipolarSpheres

        private
        
        ! Particles
        
        real(DP), dimension(:, :), allocatable, public :: M !< moments of all particles
        
        ! Monte-Carlo
        integer :: structure_iStep
        real(DP) :: dm !< rotation
        real(DP) :: dmSave
        real(DP) :: dmMax
        real(DP) :: rejRotFix
        integer :: NadaptRot

        ! Potential
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :), allocatable :: Epot_real_tab !< tabulation : real short-range
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: Epot_reci_weight
        integer :: NwaveVectors
        complex(DP), dimension(Dim, -Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: structure
        complex(DP), dimension(:, :), allocatable :: potential
        real(C_double), dimension(Dim) :: moduli_drifted
        real(C_double), dimension(Dim) :: moduli_nfft
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => DipolarSpheres_construct
        procedure :: destroy => DipolarSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: report => DipolarSpheres_report
        
        !> Take a snap shot of the configuration : orientations
        procedure :: snapShot_M => DipolarSpheres_snapShot_M
        procedure :: C_snapShot => DipolarSpheres_C_snapShot
        
        !> Adapt the displacement dx during thermalisation
        procedure :: adaptDm => DipolarSpheres_adaptDm
        procedure :: definiteDm => DipolarSpheres_definiteDm
        procedure :: getDm => DipolarSpheres_getDm
        procedure :: getNadaptRot => DipolarSpheres_getNadaptRot
        
        procedure :: getStructure_iStep => DipolarSpheres_getStructure_iStep
        
        !> Potential energy
        !>     Real
        procedure :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure :: Epot_real_interpol => DipolarSpheres_Epot_real_interpol
        procedure :: Epot_real_pair => DipolarSpheres_Epot_real_pair
        procedure :: Epot_real => DipolarSpheres_Epot_real
        !>     Reciprocal
        procedure :: Epot_reci_init => DipolarSpheres_Epot_reci_init
        procedure :: Epot_reci_weight_init => DipolarSpheres_Epot_reci_weight_init
        procedure :: Epot_reci_structure_init => DipolarSpheres_Epot_reci_structure_init
        procedure :: Epot_reci_structure_reInit => DipolarSpheres_Epot_reci_structure_reInit
        procedure :: Epot_reci_potential_init => DipolarSpheres_Epot_reci_potential_init
        procedure :: Epot_reci => DipolarSpheres_Epot_reci
        !>     Self
        procedure :: Epot_self_solo => DipolarSpheres_Epot_self_solo
        procedure :: Epot_self => DipolarSpheres_Epot_self
        !>     (Other)
        procedure :: Epot_neigh => DipolarSpheres_Epot_neigh
        procedure :: Epot_conf => DipolarSpheres_Epot_conf
        procedure :: consistTest => DipolarSpheres_consistTest
        
        !> Monte-Carlo
        procedure :: move => DipolarSpheres_move
        procedure :: rotate => DipolarSpheres_rotate
        procedure :: widom => DipolarSpheres_widom
        
    end type DipolarSpheres
    
contains

    subroutine DipolarSpheres_construct(this, shared_cell_Lsize, shared_rCut)
    
        class(DipolarSpheres), intent(out) :: this
        real(DP), dimension(:), intent(in) :: shared_cell_Lsize
        real(DP), intent(in) :: shared_rCut
        
        this%name = "dipol"
    
        ! Particles
        this%radius = dipol_radius
        this%rMin = dipol_rMin
        this%Ncol = dipol_Ncol
        allocate(this%X(Dim, this%Ncol))
        allocate(this%M(Dim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = dipol_snap_factor
        
        ! Monte-Carlo
        this%structure_iStep = dipol_structure_iStep
        this%dx = dipol_dx
        this%dxSave = this%dx
        this%rejFix = dipol_rejFix
        this%Nadapt = dipol_Nadapt
        
        this%dm = dipol_dm
        this%dmSave = this%dm
        this%dmMax = dipol_dmMax
        this%rejRotFix = dipol_rejRotFix
        this%NadaptRot = dipol_NadaptRot
        
        this%Nwidom = dipol_Nwidom
        
        ! Potential
        this%rCut = dipol_rCut
        this%dr = dipol_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%alpha = dipol_alpha        
        allocate(this%Epot_real_tab(this%iMin:this%iCut, 2))
        call this%Epot_real_init()

        allocate(this%potential(Dim, this%Ncol))
        call this%Epot_reci_weight_init()
        call C_Epot_reci_nfft_init(int(this%Ncol, C_int))
        
        ! Neighbours : same kind
        call this%same%construct(dipol_cell_Lsize, this%rCut)
        call this%same%alloc_cells()
        call this%same%ini_cell_neighs()
        ! Neighbours : other kind
        call this%mix%construct(shared_cell_Lsize, shared_rCut)
        call this%mix%alloc_cells()
        call this%mix%ini_cell_neighs()
    
    end subroutine DipolarSpheres_construct
    
    subroutine DipolarSpheres_destroy(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        if (allocated(this%X)) then
            deallocate(this%X)
        end if
        
        if (allocated(this%M)) then
            deallocate(this%M)
        end if
        
        if (allocated(this%Epot_real_tab)) then
            deallocate(this%Epot_real_tab)
        endif

        if (allocated(this%potential))then
            deallocate(this%potential)
        end if
        call C_Epot_reci_nfft_finalize()
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine DipolarSpheres_destroy
    
    !> Report
    
    subroutine DipolarSpheres_report(this, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt
        write(report_unit, *) "    Structure_iStep = ", this%structure_iStep

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        
        write(report_unit, *) "    same_cell_coordMax(:) = ", this%same%cell_coordMax(:)
        write(report_unit, *) "    same_cell_Lsize(:) = ", this%same%cell_Lsize(:)
        write(report_unit, *) "    mix_cell_coordMax(:) = ", this%mix%cell_coordMax(:)
        write(report_unit, *) "    mix_cell_Lsize(:) = ", this%mix%cell_Lsize(:)
        
    end subroutine DipolarSpheres_report
    
    !> Configuration state : orientations
      
    subroutine DipolarSpheres_snapShot_M(this, iStep, snap_unit)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%M(:, iCol)
            end do
            
        end if

    end subroutine DipolarSpheres_snapShot_M
    
    subroutine DipolarSpheres_C_snapShot(this)
        
        class(DipolarSpheres), intent(in) :: this
        
        call C_snapShot(int(this%Ncol, C_int))

    end subroutine DipolarSpheres_C_snapShot
    
    !> Adaptation of dm during the thermalisation
    
    subroutine DipolarSpheres_adaptDm(this, rej)
    
        class(DipolarSpheres), intent(inout) :: this
        real(DP), intent(in) :: rej
        
        real(DP), parameter :: eps_dm = 0.05_DP
        real(DP), parameter :: eps_rej = 0.1_DP * eps_dm
        real(DP), parameter :: more = 1._DP+eps_dm
        real(DP), parameter :: less = 1._DP-eps_dm

        
        if (rej < this%rejRotFix - eps_rej) then
        
            this%dm = this%dm * more
  
            if (this%dm > this%dmMax) then
                this%dm = this%dmMax
            end if
            
        else if (rej > this%rejRotFix + eps_rej) then
        
            this%dm = this%dm * less
            
        end if
    
    end subroutine DipolarSpheres_adaptDm
    
    subroutine DipolarSpheres_definiteDm(this, rej, report_unit)
    
        class(DipolarSpheres), intent(inout) :: this    
        real(DP), intent(in) :: rej
        integer, intent(in) :: report_unit
        
        if (rej == 0._DP) then
            write(error_unit, *) this%name, " :    Warning : dm adaptation problem."
            this%dm = this%dmSave
            write(error_unit, *) "default dm :", this%dm
        end if
        
        if (this%dm > this%dmMax) then
            write(error_unit, *) this%name, " :   Warning : dm too big."
            this%dm = this%dmMax
            write(error_unit, *) "big dm :", this%dm
        end if
        
        write(output_unit, *) this%name, " :    Thermalisation : over (rotation)"
        
        write(report_unit, *) "Rotation :"
        write(report_unit, *) "    dm = ", this%dm
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(rej-this%rejRotFix)/this%rejRotFix
    
    end subroutine DipolarSpheres_definiteDm
    
    !> Accessor : dm
    
    function DipolarSpheres_getDm(this) result(getDm)
        
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: getDm
        
        getDm = this%dm
        
    end function DipolarSpheres_getDm
    
    !> Accessor : Nadapt
    
    function DipolarSpheres_getNadaptRot(this) result(getNadaptRot)
    
        class(DipolarSpheres), intent(in) :: this        
        integer :: getNadaptRot
        
        getNadaptRot = this%NadaptRot
        
    end function DipolarSpheres_getNadaptRot
    
    !> Accessor : structure_iStep
    
    function DipolarSpheres_getStructure_iStep(this) result (getStructure_iStep)
    
        class(DipolarSpheres), intent(in) :: this
        integer :: getStructure_iStep
    
        getStructure_iStep = this%structure_iStep
        
    end function DipolarSpheres_getStructure_iStep

    ! Real -----------------------------------------------------------------------------------------
    
    !> Potential energy : real part
    !> Initialisation
    !> \f[ B(r) = \frac{\mathrm{erfc}(\alpha r)}{r^3} + 
    !>           2\frac{\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 r^2}}{r^2} \f]
    !> \f[ C(r) = 3\frac{\mathrm{erfc}(\alpha r)}{r^5} +
    !>            2\frac{\alpha}{\sqrt{\pi}}(2\alpha^2 + \frac{3}{r^2})
    !>                                     \frac{e^{-\alpha^2 r^2}}{r^2} \f]
    
    subroutine DipolarSpheres_Epot_real_init(this)
    
        class(DipolarSpheres), intent(inout) :: this

        integer :: i
        real(DP) :: r_i
        real(DP) :: alpha
        
        alpha = this%alpha
       
        ! cut
        do i = this%iMin, this%iCut
        
            r_i = real(i, DP)*this%dr
            
            this%Epot_real_tab(i, 1) = erfc(alpha*r_i)/r_i**3 + &
                                  2._DP*alpha/sqrt(PI) * exp(-alpha**2*r_i**2) / r_i**2
                                 
            this%Epot_real_tab(i, 2) = 3._DP*erfc(alpha*r_i)/r_i**5 + &
                                  2._DP*alpha/sqrt(PI) * (2_DP*alpha**2+3._DP/r_i**2) * &
                                                         exp(-alpha**2*r_i**2) / r_i**2
                                    
        end do
        
        ! shift        
        this%Epot_real_tab(:, 1) = this%Epot_real_tab(:, 1) - this%Epot_real_tab(this%iCut, 1)
        this%Epot_real_tab(:, 2) = this%Epot_real_tab(:, 2) - this%Epot_real_tab(this%iCut, 2)

    end subroutine DipolarSpheres_Epot_real_init
    
    !> Potential energy : real part
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
    
    !> Potential energy : real part
    !> Linear interpolation

    function DipolarSpheres_Epot_real_interpol(this, r) result(Epot_real_interpol)
        
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
    
    !> Potential energy : real part
    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B(r_{ij}) - 
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C(r_{ij}) \f]
    
    function DipolarSpheres_Epot_real_pair(this, mCol_i, mCol_j, rVec_ij, r_ij) &
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
    
    !> Total real energy
    
    function DipolarSpheres_Epot_real(this) result(Epot_real)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_real
        
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij
        integer :: iCol, jCol
    
        Epot_real = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    rVec_ij = distVec(this%X(:, iCol), this%X(:, jCol))
                    r_ij = norm2(rVec_ij)
                    
                    Epot_real = Epot_real + &
                                this%Epot_real_pair(this%M(:, iCol), this%M(:, jCol), rVec_ij, r_ij)
                    
                end if
            end do
        end do
        
        Epot_real = 0.5_DP*Epot_real
    
    end function DipolarSpheres_Epot_real

    ! Reciprocal -----------------------------------------------------------------------------------

    subroutine DipolarSpheres_Epot_reci_init(this)

        class(DipolarSpheres), intent(inout) :: this

        call this%Epot_reci_structure_init()
        call this%Epot_reci_potential_init()

    end subroutine
    
    !> \f[ f(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}
    !>      {\alpha^2} \sum_d \frac{k_d^2}{L_d}}}{\sum_d \frac{k_d^2}{L_d}} \f]
    
    subroutine DipolarSpheres_Epot_reci_weight_init(this)
        
        class(DipolarSpheres), intent(inout) :: this
        
        integer :: kx, ky, kz
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: kOverL

        this%NwaveVectors = 0

        do kz = -Kmax(3), Kmax(3)
            waveVector(3) = real(kz, DP)
        
        do ky = -Kmax(2), Kmax(2)
            waveVector(2) = real(ky, DP)
        
        do kx = -Kmax(1), Kmax(1)
            waveVector(1) = real(kx, DP)

            if (norm2(waveVector) /= 0) then
            
                kOverL = norm2(waveVector(:)/Lsize(:))

                this%Epot_reci_weight(kx, ky, kz) = exp(-PI**2/this%alpha**2 * kOverL**2) / kOverL**2

                this%NwaveVectors = this%NwaveVectors + 1

            else

                this%Epot_reci_weight(kx, ky, kz) = 0._DP

            end if

        end do
            
        end do
        
        end do
        
    end subroutine DipolarSpheres_Epot_reci_weight_init

    subroutine DipolarSpheres_Epot_reci_structure_init(this)

        class(DipolarSpheres), intent(inout) :: this

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3

        real(DP), dimension(Dim) :: xColOverL, mColOverL
        integer :: kx, ky, kz
        integer :: iCol

        this%structure(:, :, :, :) = cmplx(0._DP, 0._DP, DP)

        do iCol = 1, this%Ncol
        
            xColOverL(:) = this%X(:, iCol)/Lsize(:)
            mColOverL(:) = this%M(:, iCol)/Lsize(:)
            
            call fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
        
            do kz = -Kmax(3), Kmax(3)
            do ky = -Kmax(2), Kmax(2)
            do kx = -Kmax(1), Kmax(1)            
            
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)
                          
                this%structure(:, kx, ky, kz) = this%structure(:, kx, ky, kz) + &
                                                cmplx(mColOverL(:), 0._DP, DP) * exp_IkxCol
            
            end do
            end do
            end do
            
        end do

        do kz = -Kmax(3), Kmax(3)
        do ky = -Kmax(2), Kmax(2)
        do kx = -Kmax(1), Kmax(1)

            write(*, *) kx, ky, kz, &
                        real(this%structure(:, kx, ky, kz)), aimag(this%structure(:, kx, ky, kz))

        end do
        end do
        end do

    end subroutine DipolarSpheres_Epot_reci_structure_init
    
    subroutine DipolarSpheres_Epot_reci_structure_reInit(this, iStep, moduli_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: moduli_unit
        
        call C_Epot_reci_structure_moduli(this%moduli_drifted)
        call C_Epot_reci_structure()
        call C_Epot_reci_structure_moduli(this%moduli_nfft)
        
        write(moduli_unit, *) iStep, abs(this%moduli_nfft(:)-this%moduli_drifted(:))
    
    end subroutine DipolarSpheres_Epot_reci_structure_reInit
    
    subroutine DipolarSpheres_Epot_reci_potential_init(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        complex(DP) :: conjg_exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3
        complex(DP) :: k_dot_structure
        
        real(DP), dimension(Dim) :: xColOverL
        real(DP), dimension(Dim) :: waveVector
        integer :: kx, ky, kz
        integer :: iCol
        
        this%potential(:, :) = cmplx(0._DP, 0._DP, DP)
        
        do iCol = 1, this%Ncol
        
            xColOverL(:) = this%X(:, iCol)/Lsize(:)
            
            call fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
        
            do kz = -Kmax(3), Kmax(3)          
              
                waveVector(3) = real(kz, DP)
            
            do ky = -Kmax(2), Kmax(2)      
                  
                waveVector(2) = real(ky, DP)
            
            do kx = -Kmax(1), Kmax(1)      
                  
                waveVector(1) = real(kx, DP)
            
                conjg_exp_IkxCol = conjg(exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz))
                
                k_dot_structure = dot_product(cmplx(waveVector, 0._DP, DP), &
                                              this%structure(:, kx, ky, kz))
                
                this%potential(:, iCol) = this%potential(:, iCol) + &
                                          cmplx(waveVector(:), 0._DP, DP) * &
                                          cmplx(this%Epot_reci_weight(kx, ky, kz), 0._DP, DP) * &
                                          k_dot_structure * conjg_exp_IkxCol
                
            end do
            end do
            end do
            
        end do
    
    end subroutine DipolarSpheres_Epot_reci_potential_init
    
    !> Total reciprocal energy
    
    function DipolarSpheres_Epot_reci(this) result(Epot_reci)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci
        
        integer :: jCol
        real(DP), dimension(Dim) :: mColOverL
        real(DP), dimension(Dim) :: real_potential

        Epot_reci = 0._DP

        do jCol = 1, this%Ncol

            mColOverL(:) = this%M(:, jCol)/Lsize(:)
            real_potential(:) = real(this%potential(:, jCol), DP)

            Epot_reci = Epot_reci + dot_product(mColOverL, real_potential)

            write(*, *) "Epot_reci :", jCol, Epot_reci
        
        end do
        
        Epot_reci = 2._DP*PI/Volume * Epot_reci
        
    end function DipolarSpheres_Epot_reci
    
    ! Self -----------------------------------------------------------------------------------------
    
    !> Self energy of 1 dipole
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    function DipolarSpheres_Epot_self_solo(this, mCol) result(Epot_self_solo)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: Epot_self_solo
        
        Epot_self_solo = 2._DP/3._DP * this%alpha**3/sqrt(PI) * dot_product(mCol, mCol)
    
    end function DipolarSpheres_Epot_self_solo
    
    !> Total self energy
    !> \f[ \frac{2}{3}\frac{\alpha^3}{\sqrt{\pi}} \sum_i \vec{\mu}_i\cdot\vec{\mu}_i \f]
    
    function DipolarSpheres_Epot_self(this) result(Epot_self)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_self

        integer :: iCol        
        
        Epot_self = 0._DP
        do iCol = 1, this%Ncol
            Epot_self = Epot_self + this%Epot_self_solo(this%M(:, iCol))
        end do
        
    end function DipolarSpheres_Epot_self
    
    ! ----------------------------------------------------------------------------------------------
    
    !> Real potential energy : short-range
    
    subroutine DipolarSpheres_Epot_neigh(this, iCol, xCol, mCol, iCell, overlap, energ)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(:), intent(in) :: xCol, mCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNeigh,  iCell_neigh
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij        
        
        type(Link), pointer :: current => null(), next => null()
        
        
        overlap = .false.
        energ = 0._DP
    
        do iNeigh = 1, cell_neighs_nb
        
            iCell_neigh = this%same%cell_neighs(iNeigh, iCell)
            current => this%same%cellsBegin(iCell_neigh)%particle%next            
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= iCol) then
                
                    rVec_ij = distVec(xCol(:), this%X(:, current%iCol))
                    r_ij = norm2(rVec_ij)
                    
                    if (r_ij < this%rMin) then
                        overlap = .true.
                        return
                    end if
                    
                    energ = energ + this%Epot_real_pair(mCol, this%M(:, current%iCol), &
                                                        rVec_ij, r_ij)
       
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine DipolarSpheres_Epot_neigh
    
    !> Particle move
    
    subroutine DipolarSpheres_move(this, iOld, other, mix, same_Epot, mix_Epot, Nrej)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: same_Epot, mix_Epot
        integer, intent(inout) :: Nrej
        
        real(DP), dimension(Dim) :: xRand
        logical :: overlap
        real(DP), dimension(Dim) :: xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: dEpot, same_dEpot, mix_dEpot
        real(DP) :: same_eNew, same_eOld
        real(DP) :: mix_eNew, mix_eOld
        real(DP) :: rand
        
        real(C_double) :: C_Epot
        real(C_double), dimension(Dim) :: C_xNew
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = this%X(:, iOld) + this%dx(:) * (xRand(:)-0.5_DP)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        mix_iCellNew = this%mix%position_to_cell(xNew)
        call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%X, overlap, mix_eNew)
            
        if (.not. overlap) then
        
            same_iCellNew = this%same%position_to_cell(xNew)
            call this%Epot_neigh(iOld, xNew, this%M(:, iOld), same_iCellNew, overlap, same_eNew)
                        
            if (.not. overlap) then
                
                ! Real
                same_iCellOld = this%same%position_to_cell(this%X(:, iOld))
                call this%Epot_neigh(iOld, this%X(:, iOld), this%M(:, iOld), same_iCellOld, &
                                     overlap, same_eOld)
                ! Reci
                C_xNew(:) = real(xNew(:)/Lsize(:), C_double) - 0.5_c_double
                C_Epot = C_Epot_reci_move(int(iOld-1, C_int), C_xNew, real(Volume, C_double))
                
                same_dEpot = (same_eNew - same_eOld) + real(C_Epot, DP)
                    
                mix_iCellOld = this%mix%position_to_cell(this%X(:, iOld))
                call mix%Epot_neigh(this%X(:, iOld), mix_iCellOld, this%mix, other%X, overlap, &
                                    mix_eOld)
                mix_dEpot = mix_eNew - mix_eOld
                
                dEpot = same_dEpot + mix_dEpot
                
                call random_number(rand)            
                if (rand < exp(-dEpot/Tstar)) then
                
                    this%X(:, iOld) = xNew(:)
                    call C_Epot_reci_updateX(int(iOld-1, C_int), C_xNew)
                    
                    same_Epot = same_Epot + same_dEpot
                    mix_Epot = mix_Epot + mix_dEpot
                    
                    if (same_iCellOld /= same_iCellNew) then
                        call this%same%remove_cell_col(iOld, same_iCellOld)
                        call this%same%add_cell_col(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then
                        call other%mix%remove_cell_col(iOld, mix_iCellOld)
                        call other%mix%add_cell_col(iOld, mix_iCellNew)
                    end if
                    
                else
                    Nrej = Nrej + 1
                end if
         
            else
                Nrej = Nrej + 1
            end if            
            
        else        
            Nrej = Nrej + 1
        end if
    
    end subroutine DipolarSpheres_move
    
    subroutine DipolarSpheres_rotate(this, iOld, Epot, Nrej)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        real(DP), intent(inout) :: Epot
        integer, intent(inout) :: Nrej
        
        real(DP) :: rand
        real(DP), dimension(Dim) :: mNew
        real(DP) :: dEpot, real_dEpot, self_dEpot
        real(DP) :: real_eNew, real_eOld
        integer :: iCell
        logical :: overlap
        
        real(C_double) :: C_Epot
        real(C_double), dimension(Dim) :: C_mNew
        
        mNew(:) = this%M(:, iOld)
        call markov_surface(mNew, this%dm)
        
        C_mNew(:) = real(mNew(:)/Lsize(:), C_double)
        C_Epot = C_Epot_reci_rotate(int(iOld-1, C_int), C_mNew, real(Volume, C_double))
        
        iCell = this%same%position_to_cell(this%X(:, iOld))
        call this%Epot_neigh(iOld, this%X(:, iOld), mNew, iCell, overlap, real_eNew)
        call this%Epot_neigh(iOld, this%X(:, iOld), this%M(:, iOld), iCell, overlap, real_eOld)        
        real_dEpot = real_eNew - real_eOld
        
        self_dEpot = this%Epot_self_solo(mNew) - this%Epot_self_solo(this%M(:, iOld))
        
        dEpot = real(C_Epot, DP) + real_dEpot - self_dEpot
        call random_number(rand)
        if (rand < exp(-dEpot/Tstar)) then
        
            this%M(:, iOld) = mNew(:)
            call C_Epot_reci_updateM(int(iOld-1, C_int), C_mNew)
            
            Epot = Epot + dEpot
            
        else
            Nrej = Nrej + 1
        end if
    
    end subroutine DipolarSpheres_rotate
    
    !> Widom's method : with other type ?

    subroutine DipolarSpheres_widom(this, other_X, mix, activ)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: other_X
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ 
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest, xRand
        real(DP), dimension(Dim) :: mTest
        integer :: same_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: enTest, same_enTest, mix_enTest
        
        real(C_double) :: C_Epot
        real(C_double), dimension(Dim) :: C_xTest, C_mTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)
            
            mix_iCellTest = this%mix%position_to_cell(xTest)
            call mix%Epot_neigh(xTest, mix_iCellTest, this%mix, other_X, overlap, mix_enTest)
            
            if (.not. overlap) then
            
                mTest(:) = random_surface()
                               
                same_iCellTest = this%same%position_to_cell(xTest)               
                call this%Epot_neigh(0, xTest, mTest, same_iCellTest, overlap, same_enTest)
                
                if (.not. overlap) then
                
                    ! Reci
                    C_xTest(:) = real(xTest(:)/Lsize(:), C_double) - 0.5_c_double                    
                    C_mTest(:) = real(mTest(:)/Lsize(:), C_double)
                    
                    C_Epot = C_Epot_reci_test(C_xTest, C_mTest, real(Volume, C_double))
                
                    enTest = same_enTest + mix_enTest + real(C_Epot, DP) - this%Epot_self_solo(mTest)
                    widTestSum = widTestSum + exp(-enTest/Tstar)
                    
                end if
            
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine DipolarSpheres_widom

    !> Total potential energy
    
    function DipolarSpheres_Epot_conf(this) result(Epot_conf)
    
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
        
        real(DP) :: Epot_real, Epot_reci, Epot_self
        
        Epot_real = this%Epot_real()
        Epot_reci = this%Epot_reci()
        Epot_self = this%Epot_self()
        
        Epot_conf = Epot_real + Epot_reci - Epot_self
        
        write(output_unit, *) "Epot_real", Epot_real
        write(output_unit, *) "Epot_reci", Epot_reci
        write(output_unit, *) "Epot_self", Epot_self
    
    end function DipolarSpheres_Epot_conf
    
    !> Consistency test 
    
    subroutine DipolarSpheres_consistTest(this, Epot, report_unit)
    
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
    
    end subroutine DipolarSpheres_consistTest

end module class_dipolarSpheres