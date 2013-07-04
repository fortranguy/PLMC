!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, consist_tiny
use data_constants, only : PI
use data_cell, only : Dim, Lsize, kMax, Volume
use data_particles, only : dipol_radius, dipol_rMin, dipol_Ncol
use data_mc, only : Temperature, dipol_structure_iStep, dipol_deltaX, dipol_rejectFix, dipol_Nadapt, &
                    dipol_deltaM, dipol_deltaMmax, dipol_rejectRotFix, dipol_NadaptRot, dipol_Nwidom
use data_potentiel, only : dipol_rCut, dipol_dr, dipol_alpha
use data_neighbours, only : cell_neighs_nb, dipol_cell_Lsize
use data_distrib, only : dipol_snap_factor
use mod_physics, only : dist, distVec, random_surface, markov_surface, fourier
use class_neighbours
use class_mixingPotential
use class_spheres

implicit none

private

    type, extends(Spheres), public :: DipolarSpheres

        private
        
        ! Particles
        
        real(DP), dimension(:, :), allocatable, public :: orientations !< dipolar orientations 
                                                                       !< of all particles
        
        ! Monte-Carlo
        integer :: structure_iStep
        real(DP) :: deltaM !< rotation
        real(DP) :: deltaMsave
        real(DP) :: deltaMmax
        real(DP) :: rejectRotFix
        integer :: NadaptRot

        ! Potential
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :), allocatable :: Epot_real_tab !< tabulation : real short-range
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: Epot_reci_weight
        integer :: NwaveVectors
        complex(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: &
            Epot_reci_kStructure
        complex(DP), dimension(:, :), allocatable :: Epot_reci_potential
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => DipolarSpheres_construct
        procedure :: destroy => DipolarSpheres_destroy
        
        !> Print a report of the component in a file
        procedure :: PrintReport => DipolarSpheres_printReport
        
        !> Take a snap shot of the configuration : orientations
        procedure :: snapShot_orientations => DipolarSpheres_snapShot_orientations
        
        !> Adapt the rotation deltaM during thermalisation
        procedure :: adaptDeltaM => DipolarSpheres_adaptDeltaM
        procedure :: definiteDeltaM => DipolarSpheres_definiteDeltaM
        procedure :: getDeltaM => DipolarSpheres_getDeltaM
        procedure :: getNadaptRot => DipolarSpheres_getNadaptRot
        
        procedure :: getStructure_iStep => DipolarSpheres_getStructure_iStep
        
        !> Potential energy
        !>     Real
        procedure :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure :: Epot_real_interpol => DipolarSpheres_Epot_real_interpol
        procedure :: Epot_real_pair => DipolarSpheres_Epot_real_pair
        procedure :: Epot_real_overlapTest => DipolarSpheres_Epot_real_overlapTest
        procedure :: Epot_real_solo => DipolarSpheres_Epot_real_solo
        procedure :: Epot_real => DipolarSpheres_Epot_real
        !>     Reciprocal : init
        procedure :: Epot_reci_init => DipolarSpheres_Epot_reci_init
        procedure :: Epot_reci_weight_init => DipolarSpheres_Epot_reci_weight_init
        procedure :: Epot_reci_structure_init => DipolarSpheres_Epot_reci_structure_init
        procedure :: Epot_reci_structure_moduli => DipolarSpheres_Epot_reci_structure_moduli
        procedure :: Epot_reci_structure_reInit => DipolarSpheres_Epot_reci_structure_reInit
        procedure :: Epot_reci_potential_init => DipolarSpheres_Epot_reci_potential_init
        procedure :: Epot_reci_countNwaveVectors => DipolarSpheres_Epot_reci_countNwaveVectors
        !>     Reciprocal : delta
        procedure :: deltaEpot_reci_move => DipolarSpheres_deltaEpot_reci_move
        procedure :: deltaEpot_reci_move_updateStructure => &
                     DipolarSpheres_deltaEpot_reci_move_updateStructure
        procedure :: deltaEpot_reci_rotate => DipolarSpheres_deltaEpot_reci_rotate
        procedure :: deltaEpot_reci_rotate_updateStructure => &
                     DipolarSpheres_deltaEpot_reci_rotate_updateStructure
        procedure :: deltaEpot_reci_test => DipolarSpheres_deltaEpot_reci_test
        !>     Reciprocal : total
        procedure :: Epot_reci => DipolarSpheres_Epot_reci
        !>     Self
        procedure :: Epot_self_solo => DipolarSpheres_Epot_self_solo
        procedure :: Epot_self => DipolarSpheres_Epot_self
        !>     Total
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
        allocate(this%positions(Dim, this%Ncol))
        allocate(this%orientations(Dim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = dipol_snap_factor
        
        ! Monte-Carlo
        this%structure_iStep = dipol_structure_iStep
        this%deltaX = dipol_deltaX
        this%deltaXsave = this%deltaX
        this%rejectFix = dipol_rejectFix
        this%Nadapt = dipol_Nadapt
        
        this%deltaM = dipol_deltaM
        this%deltaMsave = this%deltaM
        this%deltaMmax = dipol_deltaMmax
        this%rejectRotFix = dipol_rejectRotFix
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

        allocate(this%Epot_reci_potential(Dim, this%Ncol))
        call this%Epot_reci_weight_init()
        
        ! Neighbours : same kind
        call this%same%construct(dipol_cell_Lsize, this%rCut)
        call this%same%alloc_cells()
        call this%same%cell_neighs_init()
        ! Neighbours : other kind
        call this%mix%construct(shared_cell_Lsize, shared_rCut)
        call this%mix%alloc_cells()
        call this%mix%cell_neighs_init()
    
    end subroutine DipolarSpheres_construct
    
    subroutine DipolarSpheres_destroy(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        if (allocated(this%orientations)) then
            deallocate(this%orientations)
        end if
        
        if (allocated(this%Epot_real_tab)) then
            deallocate(this%Epot_real_tab)
        endif

        if (allocated(this%Epot_reci_potential))then
            deallocate(this%Epot_reci_potential)
        end if
        
        call this%same%destroy()
        call this%mix%destroy()
    
    end subroutine DipolarSpheres_destroy
    
    !> Report
    
    subroutine DipolarSpheres_printReport(this, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    Nadapt = ", this%Nadapt

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        write(report_unit, *) "    Structure_iStep = ", this%structure_iStep
        write(report_unit, *) "    NwaveVectors = ", this%NwaveVectors
        
        write(report_unit, *) "    same_cell_coordMax(:) = ", this%same%cell_coordMax(:)
        write(report_unit, *) "    same_cell_Lsize(:) = ", this%same%cell_Lsize(:)
        write(report_unit, *) "    mix_cell_coordMax(:) = ", this%mix%cell_coordMax(:)
        write(report_unit, *) "    mix_cell_Lsize(:) = ", this%mix%cell_Lsize(:)
        
    end subroutine DipolarSpheres_printReport
    
    !> Configuration state : orientations
      
    subroutine DipolarSpheres_snapShot_orientations(this, iStep, snap_unit)
        
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%orientations(:, iCol)
            end do
            
        end if

    end subroutine DipolarSpheres_snapShot_orientations
    
    !> Adaptation of deltaM during the thermalisation
    
    subroutine DipolarSpheres_adaptDeltaM(this, reject)
    
        class(DipolarSpheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: eps_deltaM = 0.05_DP
        real(DP), parameter :: eps_reject = 0.1_DP * eps_deltaM
        real(DP), parameter :: more = 1._DP+eps_deltaM
        real(DP), parameter :: less = 1._DP-eps_deltaM

        
        if (reject < this%rejectRotFix - eps_reject) then
        
            this%deltaM = this%deltaM * more
  
            if (this%deltaM > this%deltaMmax) then
                this%deltaM = this%deltaMmax
            end if
            
        else if (reject > this%rejectRotFix + eps_reject) then
        
            this%deltaM = this%deltaM * less
            
        end if
    
    end subroutine DipolarSpheres_adaptDeltaM
    
    subroutine DipolarSpheres_definiteDeltaM(this, reject, report_unit)
    
        class(DipolarSpheres), intent(inout) :: this    
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        if (reject == 0._DP) then
            write(error_unit, *) this%name, " :    Warning : deltaM adaptation problem."
            this%deltaM = this%deltaMsave
            write(error_unit, *) "default deltaM :", this%deltaM
        end if
        
        if (this%deltaM > this%deltaMmax) then
            write(error_unit, *) this%name, " :   Warning : deltaM too big."
            this%deltaM = this%deltaMmax
            write(error_unit, *) "big deltaM :", this%deltaM
        end if
        
        write(output_unit, *) this%name, " :    Thermalisation : over (rotation)"
        
        write(report_unit, *) "Rotation :"
        write(report_unit, *) "    deltaM = ", this%deltaM
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%rejectRotFix)/this%rejectRotFix
    
    end subroutine DipolarSpheres_definiteDeltaM
    
    !> Accessor : deltaM
    
    pure function DipolarSpheres_getDeltaM(this) result(getDeltaM)
        
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: getDeltaM
        
        getDeltaM = this%deltaM
        
    end function DipolarSpheres_getDeltaM
    
    !> Accessor : Nadapt
    
    pure function DipolarSpheres_getNadaptRot(this) result(getNadaptRot)
    
        class(DipolarSpheres), intent(in) :: this        
        integer :: getNadaptRot
        
        getNadaptRot = this%NadaptRot
        
    end function DipolarSpheres_getNadaptRot
    
    !> Accessor : structure_iStep
    
    pure function DipolarSpheres_getStructure_iStep(this) result (getStructure_iStep)
    
        class(DipolarSpheres), intent(in) :: this
        integer :: getStructure_iStep
    
        getStructure_iStep = this%structure_iStep
        
    end function DipolarSpheres_getStructure_iStep

    ! Real : short-range interaction ---------------------------------------------------------------
    
    !> Initialisation : look-up (tabulation) table
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
 
    !> Overlap test of 1 particle

    subroutine DipolarSpheres_Epot_real_overlapTest(this, iCol, xCol, iCell, overlap)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol, iCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap

        integer :: iNeigh,  iCell_neigh, jCol
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij

        type(Link), pointer :: current => null(), next => null()

        overlap = .false.        

        do iNeigh = 1, cell_neighs_nb

            iCell_neigh = this%same%cell_neighs(iNeigh, iCell)
            current => this%same%cellsBegin(iCell_neigh)%particle%next
            if (.not. associated(current%next)) cycle

            do

                next => current%next

                if (current%iCol /= iCol) then

                    r_ij = dist(xCol(:), this%positions(:, current%iCol))

                    if (r_ij < this%rMin) then
                        overlap = .true.
                        return
                    end if

                end if

                if (.not. associated(next%next)) exit

                current => next

            end do

        end do

    end subroutine DipolarSpheres_Epot_real_overlapTest
    
    !> Energy of 1 dipole with others
    
    pure function DipolarSpheres_Epot_real_solo(this, iCol, xCol, mCol) result(Epot_real_solo)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol
        real(DP), dimension(:), intent(in) :: xCol, mCol
        real(DP) :: Epot_real_solo

        integer :: jCol
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij

        Epot_real_solo = 0._DP

        do jCol = 1, this%Ncol

            if (jCol /= iCol) then

                rVec_ij = distVec(xCol(:), this%positions(:, jCol))
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
        real(DP), dimension(Dim) :: rVec_ij
        real(DP) :: r_ij
    
        Epot_real = 0._DP
        
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    rVec_ij = distVec(this%positions(:, iCol), this%positions(:, jCol))
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
    
    !> Initialisation of the ``potential'' and the ``structure factor''

    subroutine DipolarSpheres_Epot_reci_init(this)

        class(DipolarSpheres), intent(inout) :: this

        call this%Epot_reci_structure_init()
        call this%Epot_reci_potential_init()

    end subroutine
    
    !> \f[ w(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}
    !>      {\alpha^2} \sum_d \frac{k_d^2}{L_d}}}{\sum_d \frac{k_d^2}{L_d}} \f]
    
    subroutine DipolarSpheres_Epot_reci_weight_init(this)
        
        class(DipolarSpheres), intent(inout) :: this
        
        integer :: kx, ky, kz
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: kOverL

        do kz = -Kmax(3), Kmax(3)
            waveVector(3) = real(kz, DP)
        
        do ky = -Kmax(2), Kmax(2)
            waveVector(2) = real(ky, DP)
        
        do kx = -Kmax(1), Kmax(1)
            waveVector(1) = real(kx, DP)

            if (norm2(waveVector) /= 0._DP) then
            
                kOverL = norm2(waveVector(:)/Lsize(:))

                this%Epot_reci_weight(kx, ky, kz) = exp(-PI**2/this%alpha**2 * kOverL**2) / kOverL**2

            else

                this%Epot_reci_weight(kx, ky, kz) = 0._DP

            end if

        end do
            
        end do
        
        end do
        
    end subroutine DipolarSpheres_Epot_reci_weight_init
    
    !> Structure factor init :
    !> \f[ \vec{S}(\vec{k}) = \sum_{i} \vec{\mu}_i e^{+i\vec{k}\cdot\vec{x}_i} \f]
    !> We will also use a restricted definition later :
    !> \f[ \vec{S}_l(\vec{k}) = \sum_{i \neq l} \vec{\mu}_i e^{+i\vec{k}\cdot\vec{x}_i} \f].

    subroutine DipolarSpheres_Epot_reci_structure_init(this)

        class(DipolarSpheres), intent(inout) :: this

        complex(DP) :: exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3

        real(DP), dimension(Dim) :: xColOverL, mColOverL
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz
        integer :: iCol

        this%Epot_reci_kStructure(:, :, :) = cmplx(0._DP, 0._DP, DP)

        do iCol = 1, this%Ncol
        
            xColOverL(:) = this%positions(:, iCol)/Lsize(:)
            mColOverL(:) = this%orientations(:, iCol)/Lsize(:)
            
            call fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
        
            do kz = -Kmax(3), Kmax(3)

                waveVector(3) = real(kz, DP)

            do ky = -Kmax(2), Kmax(2)

                waveVector(2) = real(ky, DP)

            do kx = -Kmax(1), Kmax(1)

                waveVector(1) = real(kx, DP)
            
            
                exp_IkxCol = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                k_dot_mCol = dot_product(waveVector, mColOverL)
                          
                this%Epot_reci_kStructure(kx, ky, kz) = this%Epot_reci_kStructure(kx, ky, kz) + &
                                                        cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine DipolarSpheres_Epot_reci_structure_init
    
    !> Symmetry : half wave vectors in do loop : kMax2
    
    pure function kMax2_sym(kz)

        integer, intent(in) :: kz
        integer :: kMax2_sym

        if (kz == 0) then
            kMax2_sym = 0
        else
            kMax2_sym = kMax(2)
        end if

    end function kMax2_sym
    
    !> Symmetry : half wave vectors in do loop : kMax1

    pure function kMax1_sym(ky, kz)

        integer, intent(in) :: ky, kz
        integer :: kMax1_sym

        if (ky == 0 .and. kz == 0) then
            kMax1_sym = 0
        else
            kMax1_sym = kMax(1)
        end if

    end function kMax1_sym
    
    !> To calculate the drift of the strucutre factor

    pure function DipolarSpheres_Epot_reci_structure_moduli(this) result(Epot_reci_structure_moduli)

        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci_structure_moduli

        integer :: kx, ky, kz

        Epot_reci_structure_moduli = 0._DP

        do kz = 0, Kmax(3)
            do ky = -kMax2_sym(kz), Kmax(2)
                do kx = -kMax1_sym(ky, kz), kMax(1)
                
                    Epot_reci_structure_moduli = Epot_reci_structure_moduli + &
                                                 abs(this%Epot_reci_kStructure(kx, ky, kz))

                end do
            end do
        end do

    end function DipolarSpheres_Epot_reci_structure_moduli
    
    !> Reinitialise the structure factor and print the drift
    
    subroutine DipolarSpheres_Epot_reci_structure_reInit(this, iStep, moduli_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: moduli_unit

        real(DP) :: moduli_drifted, moduli_reInit
        
        moduli_drifted = this%Epot_reci_structure_moduli()
        call this%Epot_reci_structure_init()
        moduli_reInit = this%Epot_reci_structure_moduli()
        
        write(moduli_unit, *) iStep, abs(moduli_reInit - moduli_drifted)
    
    end subroutine DipolarSpheres_Epot_reci_structure_reInit
    
    !> Potential initialisation :
    !> \f[
    !>      \vec{\phi}(\vec{x}_j) = \sum_{\vec{k}\neq\vec{0}} \vec{k} w(\alpha, \vec{k})
    !>                              (\vec{k}\cdot\vec{S}(\vec{k})) e^{-i\vec{k}\cdot\vec{x}_j}
    !> \f]
    
    subroutine DipolarSpheres_Epot_reci_potential_init(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        complex(DP) :: conjg_exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3
        complex(DP) :: k_dot_structure, Epot_reci_weight

        complex(DP), dimension(Dim) :: potential_k
        real(DP), dimension(Dim) :: xColOverL
        real(DP), dimension(Dim) :: waveVector
        integer :: kx, ky, kz
        integer :: iCol
        
        this%Epot_reci_potential(:, :) = cmplx(0._DP, 0._DP, DP)
        
        do iCol = 1, this%Ncol
        
            xColOverL(:) = this%positions(:, iCol)/Lsize(:)
            
            call fourier(xColOverL, exp_Ikx_1, exp_Ikx_2, exp_Ikx_3)
        
            do kz = -Kmax(3), Kmax(3)
              
                waveVector(3) = real(kz, DP)
            
            do ky = -Kmax(2), Kmax(2)
                  
                waveVector(2) = real(ky, DP)
            
            do kx = -Kmax(1), Kmax(1)
                  
                waveVector(1) = real(kx, DP)
            
                conjg_exp_IkxCol = conjg(exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz))
                
                k_dot_structure = this%Epot_reci_kStructure(kx, ky, kz)

                Epot_reci_weight = cmplx(this%Epot_reci_weight(kx, ky, kz), 0._DP, DP)

                potential_k(:) = cmplx(waveVector(:), 0._DP, DP) * Epot_reci_weight * k_dot_structure
                
                this%Epot_reci_potential(:, iCol) = this%Epot_reci_potential(:, iCol) + &
                                                    potential_k(:) * conjg_exp_IkxCol

                
            end do
            end do
            end do
            
        end do
        
    end subroutine DipolarSpheres_Epot_reci_potential_init

    ! Count the number of wave vectors

    subroutine DipolarSpheres_Epot_reci_countNwaveVectors(this, waveVectors_unit)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: waveVectors_unit
        
        real(DP), dimension(Dim) :: waveVector
        integer :: kx, ky, kz

        this%NwaveVectors = 0

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -kMax1_sym(ky, kz), kMax(1)

                    waveVector(1) = real(kx, DP)

                    if (norm2(waveVector) /= 0._DP) then

                        write(waveVectors_unit, *) kx, ky, kz
                        write(waveVectors_unit, *)
                        write(waveVectors_unit, *)

                        this%NwaveVectors = this%NwaveVectors + 1

                    end if

                end do

            end do

        end do

    end subroutine DipolarSpheres_Epot_reci_countNwaveVectors

    !> Move

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta M^2
    !> w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta M^2 = 2\Re[
    !>                  (\vec{\mu}_l\cdot\vec{k})
    !>                  (e^{-i\vec{k}\cdot\vec{x}^\prime_l} - e^{-i\vec{k}\cdot\vec{x}_l})
    !>                  (\vec{k}\cdot\vec{S}_l)
    !>               ]
    !> \f]

    !> Implementation :
    !> \f[
    !>  \Delta M^2 = 2(\vec{\mu_l}\cdot\vec{k})
    !>              [ \cos(\vec{k}\cdot\vec{x}^\prime_l) - \cos(\vec{k}\cdot\vec{x}_l)]
    !>              [\Re{(\vec{k}\cdot\vec{S})} - (\vec{k}\cdot\vec{\mu}_l) \cos(\vec{k}\cdot
    !>                  \vec{x}_l)] -
    !>              [-\sin(\vec{k}\cdot\vec{x}^\prime_l) + \sin(\vec{k}\cdot\vec{x}_l)]
    !>              [\Im{(\vec{k}\cdot\vec{S})} - (\vec{k}\cdot\vec{\mu}_l) \sin(\vec{k}\cdot
    !>                  \vec{x}_l)]
    !> \f]
    !>

    pure function DipolarSpheres_deltaEpot_reci_move(this, lCol, xNew) result(deltaEpot_reci_move)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: lCol
        real(DP), dimension(Dim), intent(in) :: xNew
        real(DP) :: deltaEpot_reci_move
        
        real(DP) :: deltaEpot_k

        real(DP), dimension(Dim) :: xNewOverL, xOldOverL
        real(DP), dimension(Dim) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew
        real(DP) :: cos_kxNew, sin_kxNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld
        real(DP) :: cos_kxOld, sin_kxOld

        real(DP) :: realPart1, realPart2
        
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_mCol
        complex(DP) :: k_dot_structure
        integer :: kx, ky, kz

        xNewOverL(:) = xNew(:)/Lsize(:)
        xOldOverL(:) = this%positions(:, lCol)/Lsize(:)
        
        call fourier(xNewOverL, exp_IkxNew_1, exp_IkxNew_2, exp_IkxNew_3)
        call fourier(xOldOverL, exp_IkxOld_1, exp_IkxOld_2, exp_IkxOld_3)

        mColOverL(:) = this%orientations(:, lCol)/Lsize(:)

        deltaEpot_reci_move = 0._DP

        do kz = 0, Kmax(3) ! symmetry : half wave vectors -> double Energy

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2) 

                waveVector(2) = real(ky, DP)
            
                do kx = -kMax1_sym(ky, kz), kMax(1)

                    waveVector(1) = real(kx, DP)

                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    k_dot_structure = this%Epot_reci_kStructure(kx, ky, kz)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    cos_kxNew = real(exp_IkxNew, DP)
                    sin_kxNew = aimag(exp_IkxNew)

                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                    cos_kxOld = real(exp_IkxOld, DP)
                    sin_kxOld = aimag(exp_IkxOld)

                    realPart1 = (cos_kxNew - cos_kxOld)
                    realPart1 = realPart1 * (real(k_dot_structure, DP) - k_dot_mCol * cos_kxOld)

                    realPart2 = (-sin_kxNew + sin_kxOld)
                    realPart2 = realPart2 * (aimag(k_dot_structure) - k_dot_mCol * sin_kxOld)
                    
                    deltaEpot_k = 2._DP*k_dot_mCol * (realPart1 - realPart2) * &
                                  this%Epot_reci_weight(kx, ky, kz)
                    deltaEpot_reci_move = deltaEpot_reci_move + deltaEpot_k
                                                      

                end do
            
            end do
        
        end do

        deltaEpot_reci_move = 4._DP*PI/Volume * deltaEpot_reci_move

    end function DipolarSpheres_deltaEpot_reci_move

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta \vec{S} = \vec{\mu}_l
    !>      (e^{+i\vec{k}\cdot\vec{x}^\prime_l} -e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    subroutine DipolarSpheres_deltaEpot_reci_move_updateStructure(this, lCol, xNew)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: lCol
        real(DP), dimension(Dim), intent(in) :: xNew
        
        real(DP), dimension(Dim) :: xNewOverL, xOldOverL
        real(DP), dimension(Dim) :: mColOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_mCol
        integer :: kx, ky, kz

        xNewOverL(:) = xNew(:)/Lsize(:)
        xOldOverL(:) = this%positions(:, lCol)/Lsize(:)
        
        call fourier(xNewOverL, exp_IkxNew_1, exp_IkxNew_2, exp_IkxNew_3)
        call fourier(xOldOverL, exp_IkxOld_1, exp_IkxOld_2, exp_IkxOld_3)

        mColOverL(:) = this%orientations(:, lCol)/Lsize(:)

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -kMax1_sym(ky, kz), kMax(1)

                    waveVector(1) = real(kx, DP)

                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                                                          
                    this%Epot_reci_kStructure(kx, ky, kz) = &
                        this%Epot_reci_kStructure(kx, ky, kz) + &
                        cmplx(k_dot_mCol, 0._DP, DP) * (exp_IkxNew - exp_IkxOld)

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_deltaEpot_reci_move_updateStructure
    
    !> Rotate

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta M^2
    !>                                       w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta M^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2\Re\{
    !>                  [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>                  e^{-i \vec{k} \cdot \vec{x}_l}
    !>                  (\vec{k} \cdot \vec{S}_l)
    !>               \}
    !> \f]
    
    !> Implementation :
    !> \f[
    !>  \Delta M^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2 [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>               \{
    !>                  \cos(\vec{k} \cdot \vec{x}_l)[\Re(\vec{k} \cdot \vec{S}) -
    !>                      (\vec{k} \cdot \vec{\mu}_l) \cos(\vec{k} \cdot \vec{x}_l)] +
    !>                  \sin(\vec{k} \cdot \vec{x}_l)[\Im(\vec{k} \cdot \vec{S}) -
    !>                      (\vec{k} \cdot \vec{\mu}_l) \sin(\vec{k} \cdot \vec{x}_l)]
    !>               \}
    !> \f]

    pure function DipolarSpheres_deltaEpot_reci_rotate(this, lCol, mNew) result(deltaEpot_reci_rotate)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: lCol
        real(DP), dimension(Dim), intent(in) :: mNew
        real(DP) :: deltaEpot_reci_rotate

        real(DP) :: deltaEpot_k

        real(DP), dimension(Dim) :: xColOverL
        real(DP), dimension(Dim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol
        real(DP) :: cos_kxCol, sin_kxCol

        real(DP) :: realPart, realPart1, realPart2
        
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_mNew, k_dot_mOld
        complex(DP) :: k_dot_structure
        integer :: kx, ky, kz

        xColOverL(:) = this%positions(:, lCol)/Lsize(:)
        
        call fourier(xColOverL, exp_IkxCol_1, exp_IkxCol_2, exp_IkxCol_3)

        mNewOverL(:) = mNew(:)/Lsize(:)
        mOldOverL(:) = this%orientations(:, lCol)/Lsize(:)

        deltaEpot_reci_rotate = 0._DP

        do kz = 0, Kmax(3) ! symmetry : half wave vectors -> double Energy

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)
            
                do kx = -kMax1_sym(ky, kz), kMax(1)

                    waveVector(1) = real(kx, DP)

                    k_dot_mNew = dot_product(waveVector, mNewOverL)

                    k_dot_mOld = dot_product(waveVector, mOldOverL)

                    k_dot_structure = this%Epot_reci_kStructure(kx, ky, kz)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)
                    cos_kxCol = real(exp_IkxCol, DP)
                    sin_kxCol = aimag(exp_IkxCol)

                    realPart1 = cos_kxCol * (real(k_dot_structure, DP)-k_dot_mOld*cos_kxCol)
                    realPart2 = sin_kxCol * (aimag(k_dot_structure)-k_dot_mOld*sin_kxCol)

                    realPart = realPart1 + realPart2

                    deltaEpot_k = k_dot_mNew**2 - k_dot_mOld**2
                    deltaEpot_k = deltaEpot_k + 2._DP*(k_dot_mNew - k_dot_mOld) * realPart
                    deltaEpot_k = deltaEpot_k * this%Epot_reci_weight(kx, ky, kz)
                    deltaEpot_reci_rotate = deltaEpot_reci_rotate + deltaEpot_k

                end do

            end do

        end do

        deltaEpot_reci_rotate = 4._DP*PI/Volume * deltaEpot_reci_rotate

    end function DipolarSpheres_deltaEpot_reci_rotate

    !> Update moment -> update the ``structure factor''
    !>  \f[
    !>      \Delta \vec{S} = (\vec{\mu}_l^\prime - \vec{\mu}_l) e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    subroutine DipolarSpheres_deltaEpot_reci_rotate_updateStructure(this, lCol, mNew)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: lCol
        real(DP), dimension(Dim), intent(in) :: mNew

        real(DP), dimension(Dim) :: xColOverL
        real(DP), dimension(Dim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol

        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_deltaMcol
        integer :: kx, ky, kz

        xColOverL(:) = this%positions(:, lCol)/Lsize(:)
        
        call fourier(xColOverL, exp_IkxCol_1, exp_IkxCol_2, exp_IkxCol_3)

        mNewOverL(:) = mNew(:)/Lsize(:)
        mOldOverL(:) = this%orientations(:, lCol)/Lsize(:)

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -kMax1_sym(ky, kz), kMax(1)

                    waveVector(1) = real(kx, DP)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    k_dot_deltaMcol = dot_product(waveVector, mNewOverL - mOldOverL)

                    this%Epot_reci_kStructure(kx, ky, kz) = &
                        this%Epot_reci_kStructure(kx, ky, kz) + &
                        cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_IkxCol

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_deltaEpot_reci_rotate_updateStructure

    ! Widom
    
    !> Difference of Energy 
    !> \f[ \Delta U^{N+1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}} 
    !>                          (\vec{k} \cdot \vec{\mu}_{N+1}) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot \vec{\mu}_{N+1}) + 
    !>                              2\Re[(\vec{k} \cdot \vec{S}) e^{-i \vec{k} \cdot \vec{x}_{N+1}}]
    !>                          \}
    !> \f]
    
    !> Implementation :
    !> \f[ \Delta U^{N+1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}}
    !>                          (\vec{k} \cdot \vec{\mu}_{N+1}) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot \vec{\mu}_{N+1}) +
    !>                              2 [\Re(\vec{k} \cdot \vec{S}) \cos(\vec{k} \cdot \vec{x}_{N+1}) +
    !>                                 \Im(\vec{k} \cdot \vec{S}) \sin(\vec{k} \cdot \vec{x}_{N+1})]
    !>                          \}
    !> \f]

    pure function DipolarSpheres_deltaEpot_reci_test(this, xTest, mTest) result(deltaEpot_reci_test)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(Dim), intent(in) :: xTest
        real(DP), dimension(Dim), intent(in) :: mTest
        real(DP) :: deltaEpot_reci_test
        
        real(DP) :: deltaEpot_k
        
        real(DP), dimension(Dim) :: xTestOverL
        real(DP), dimension(Dim) :: mTestOverL
        
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxTest_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxTest_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxTest_3
        complex(DP) :: exp_IkxTest
        real(DP) :: cos_kxTest, sin_kxTest
        
        real(DP) :: realPart
        
        real(DP), dimension(Dim) :: waveVector
        real(DP) :: k_dot_mTest
        complex(DP) :: k_dot_structure
        integer :: kx, ky, kz
        
        xTestOverL(:) = xTest(:)/Lsize(:)
        
        call fourier(xTestOverL, exp_IkxTest_1, exp_IkxTest_2, exp_IkxTest_3)
        
        mTestOverL(:) = mTest(:)/Lsize(:)
        
        deltaEpot_reci_test = 0._DP
        
        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -kMax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)
            
                do kx = -kMax1_sym(ky, kz), kMax(1)
                
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_mTest = dot_product(waveVector, mTestOverL)
                    
                    k_dot_structure = this%Epot_reci_kStructure(kx, ky, kz)
                                                  
                    exp_IkxTest = exp_IkxTest_1(kx) * exp_IkxTest_2(ky) * exp_IkxTest_3(kz)
                    cos_kxTest = real(exp_IkxTest, DP)
                    sin_kxTest = aimag(exp_IkxTest)
                    
                    realPart = real(k_dot_structure, DP) * cos_kxTest
                    realPart = realPart + aimag(k_dot_structure) * sin_kxTest
                    
                    deltaEpot_k = k_dot_mTest * (k_dot_mTest + 2._DP * realPart)
                    deltaEpot_k = deltaEpot_k * this%Epot_reci_weight(kx, ky, kz)
                    deltaEpot_reci_test = deltaEpot_reci_test + deltaEpot_k
                   
                end do
            
            end do
        
        end do
        
        deltaEpot_reci_test = 4._DP*PI/Volume * deltaEpot_reci_test

    end function DipolarSpheres_deltaEpot_reci_test
    
    !> Total reciprocal energy
    
    pure function DipolarSpheres_Epot_reci(this) result(Epot_reci)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci
        
        integer :: jCol
        real(DP), dimension(Dim) :: mColOverL
        real(DP), dimension(Dim) :: real_potential

        Epot_reci = 0._DP

        do jCol = 1, this%Ncol

            mColOverL(:) = this%orientations(:, jCol)/Lsize(:)
            real_potential(:) = real(this%Epot_reci_potential(:, jCol), DP)

            Epot_reci = Epot_reci + dot_product(mColOverL, real_potential)
        
        end do
        
        Epot_reci = 2._DP*PI/Volume * Epot_reci
        
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

    !> Particle move
    
    subroutine DipolarSpheres_move(this, iOld, other, mix, same_Epot, mix_Epot, Nreject)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        class(Spheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: same_Epot, mix_Epot
        integer, intent(inout) :: Nreject
        
        real(DP), dimension(Dim) :: xRand
        logical :: overlap
        real(DP), dimension(Dim) :: xNew
        integer :: same_iCellOld, same_iCellNew
        integer :: mix_iCellOld, mix_iCellNew
        real(DP) :: deltaEpot, same_deltaEpot, mix_deltaEpot
        real(DP) :: same_eNew_real, same_eOld_real
        real(DP) :: mix_eNew, mix_eOld
        real(DP) :: rand
        
        ! Random new position
        call random_number(xRand)
        xNew(:) = this%positions(:, iOld) + this%deltaX(:) * (xRand(:)-0.5_DP)
        xNew(:) = modulo(xNew(:), Lsize(:))
        
        mix_iCellNew = this%mix%position_to_cell(xNew)
        call mix%Epot_neigh(xNew, mix_iCellNew, this%mix, other%positions, overlap, mix_eNew)
            
        if (.not. overlap) then
        
            same_iCellNew = this%same%position_to_cell(xNew)
            call this%Epot_real_overlapTest(iOld, xNew, same_iCellNew, overlap)
                        
            if (.not. overlap) then
                
                ! Real
                same_iCellOld = this%same%position_to_cell(this%positions(:, iOld))
                same_eNew_real = this%Epot_real_solo(iOld, xNew, this%orientations(:, iOld))
                same_eOld_real = this%Epot_real_solo(iOld, this%positions(:, iOld), &
                                                     this%orientations(:, iOld))
                
                same_deltaEpot = (same_eNew_real-same_eOld_real) + this%deltaEpot_reci_move(iOld, xNew)
                    
                mix_iCellOld = this%mix%position_to_cell(this%positions(:, iOld))
                call mix%Epot_neigh(this%positions(:, iOld), mix_iCellOld, this%mix, other%positions, &
                                    overlap, mix_eOld)
                mix_deltaEpot = mix_eNew - mix_eOld
                
                deltaEpot = same_deltaEpot + mix_deltaEpot
                
                call random_number(rand)            
                if (rand < exp(-deltaEpot/Temperature)) then

                    call this%deltaEpot_reci_move_updateStructure(iOld, xNew)
                    this%positions(:, iOld) = xNew(:)
                    
                    same_Epot = same_Epot + same_deltaEpot
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (same_iCellOld /= same_iCellNew) then
                        call this%same%remove_col_from_cell(iOld, same_iCellOld)
                        call this%same%add_col_to_cell(iOld, same_iCellNew)
                    end if
                    
                    if (mix_iCellOld /= mix_iCellNew) then
                        call other%mix%remove_col_from_cell(iOld, mix_iCellOld)
                        call other%mix%add_col_to_cell(iOld, mix_iCellNew)
                    end if
                    
                else
                    Nreject = Nreject + 1
                end if
         
            else
                Nreject = Nreject + 1
            end if            
            
        else        
            Nreject = Nreject + 1
        end if
    
    end subroutine DipolarSpheres_move
    
    subroutine DipolarSpheres_rotate(this, iOld, Epot, Nreject)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iOld
        real(DP), intent(inout) :: Epot
        integer, intent(inout) :: Nreject
        
        real(DP) :: rand
        real(DP), dimension(Dim) :: mNew
        real(DP) :: deltaEpot, deltaEpot_real, deltaEpot_self
        real(DP) :: real_eNew, real_eOld
        integer :: iCell
        logical :: overlap
        
        mNew(:) = this%orientations(:, iOld)
        call markov_surface(mNew, this%deltaM)
        
        iCell = this%same%position_to_cell(this%positions(:, iOld))
        real_eNew = this%Epot_real_solo(iOld, this%positions(:, iOld), mNew)
        real_eOld = this%Epot_real_solo(iOld, this%positions(:, iOld), this%orientations(:, iOld))
        deltaEpot_real = real_eNew - real_eOld        
        
        deltaEpot_self = this%Epot_self_solo(mNew) - this%Epot_self_solo(this%orientations(:, iOld))
        
        deltaEpot = deltaEpot_real + this%deltaEpot_reci_rotate(iOld, mNew) - deltaEpot_self
        
        call random_number(rand)
        if (rand < exp(-deltaEpot/Temperature)) then
        
            call this%deltaEpot_reci_rotate_updateStructure(iOld, mNew)
            this%orientations(:, iOld) = mNew(:)
            
            Epot = Epot + deltaEpot
            
        else
            Nreject = Nreject + 1
        end if
    
    end subroutine DipolarSpheres_rotate
    
    !> Widom's method

    subroutine DipolarSpheres_widom(this, other_positions, mix, activ)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:, :), intent(in) :: other_positions
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ 
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Dim) :: xTest, xRand
        real(DP), dimension(Dim) :: mTest
        integer :: same_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: enTest
        real(DP) :: same_enTest, same_enTest_real
        real(DP) :: mix_enTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)
            
            mix_iCellTest = this%mix%position_to_cell(xTest)
            call mix%Epot_neigh(xTest, mix_iCellTest, this%mix, other_positions, overlap, mix_enTest)
            
            if (.not. overlap) then
                               
                same_iCellTest = this%same%position_to_cell(xTest)
                call this%Epot_real_overlapTest(0, xTest, same_iCellTest, overlap)
                
                if (.not. overlap) then
                
                    mTest(:) = random_surface()
                    same_enTest_real = this%Epot_real_solo(0, xTest, mTest)
                                        
                    same_enTest = same_enTest_real + this%deltaEpot_reci_test(xTest, mTest) - &
                                  this%Epot_self_solo(mTest)
                
                    enTest = same_enTest + mix_enTest
                    widTestSum = widTestSum + exp(-enTest/Temperature)
                    
                end if
            
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine DipolarSpheres_widom

    !> Total potential energy
    
    pure function DipolarSpheres_Epot_conf(this) result(Epot_conf)
    
        class(DipolarSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
        
        Epot_conf = this%Epot_real() + this%Epot_reci() - this%Epot_self()
    
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

