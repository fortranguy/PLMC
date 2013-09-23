!> \brief Description of the DipolarSpheres class

module class_dipolarSpheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, consist_tiny
use data_constants, only : PI
use data_box, only : Ndim, Lsize, Kmax, Volume, out_permittivity
use data_particles, only : dipol_Ncol
use data_monteCarlo, only : Temperature, dipol_move_delta, dipol_move_rejectFix, dipol_move_Nadapt, &
                            dipol_rotate_delta, dipol_rotate_deltaMax, dipol_rotate_rejectFix, &
                            dipol_rotate_Nadapt, dipol_Nwidom, dipol_structure_iStep, &
                            dipol_totalMoment_iStep
use data_potential, only : dipol_rCut, dipol_dr, dipol_alpha
use data_neighbourCells, only : NnearCell, dipol_cell_size
use data_distribution, only : dipol_snap_factor
use module_physics, only : distVec_PBC, dist_PBC, random_surface, markov_surface, Kmax1_sym, &
                           Kmax2_sym, fourier
use class_observables
use class_neighbourCells
use class_mixingPotential
use class_hardSpheres

implicit none

private

    type, extends(HardSpheres), public :: DipolarSpheres

        private
        
        ! Particles
        
        real(DP), dimension(:, :), allocatable, public :: orientations !< dipolar orientations 
                                                                       !< of all particles
        
        ! Monte-Carlo
        real(DP) :: rotate_delta !< rotation
        real(DP) :: rotate_deltaSave
        real(DP) :: rotate_deltaMax
        real(DP) :: rotate_rejectFix
        integer :: rotate_Nadapt
        integer :: structure_iStep
        integer :: totalMoment_iStep

        ! Potential
        real(DP)  :: dr !< discretisation step
        integer :: iMin !< minimum index of tabulation : minimum distance
        integer :: iCut !< maximum index of tabulation : until potential cut
        real(DP) :: alpha !< coefficient of Ewald summation
        real(DP), dimension(:, :), allocatable :: Epot_real_tab !< tabulation : real short-range
        real(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: Epot_reci_weight
        integer :: NwaveVectors
        complex(DP), dimension(-Kmax(1):Kmax(1), -Kmax(2):Kmax(2), -Kmax(3):Kmax(3)) :: structureFactor
        complex(DP), dimension(:, :), allocatable :: Epot_reci_potential
        real(DP), dimension(Ndim) :: totalMoment
        
    contains

        !> Construction and destruction of the class
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
        procedure :: get_rotate_Nadapt => DipolarSpheres_get_rotate_Nadapt
        
        procedure :: get_structure_iStep => DipolarSpheres_get_structure_iStep
        procedure :: get_totalMoment_iStep => DipolarSpheres_get_totalMoment_iStep
        
        !> Potential energy
        !>     Real
        procedure, private :: Epot_real_init => DipolarSpheres_Epot_real_init
        procedure :: Epot_real_print => DipolarSpheres_Epot_real_print
        procedure, private :: Epot_real_interpol => DipolarSpheres_Epot_real_interpol
        procedure :: Epot_real_pair => DipolarSpheres_Epot_real_pair
        procedure :: Epot_real_test_overlap => DipolarSpheres_Epot_real_test_overlap
        procedure, private :: Epot_real_solo => DipolarSpheres_Epot_real_solo
        procedure, private :: Epot_real => DipolarSpheres_Epot_real
        !>     Reciprocal : init
        procedure :: Epot_reci_init => DipolarSpheres_Epot_reci_init
        procedure, private :: Epot_reci_init_weight => DipolarSpheres_Epot_reci_init_weight
        procedure, private :: Epot_reci_init_structure => DipolarSpheres_Epot_reci_init_structure
        procedure, private :: Epot_reci_get_structure_modulus => &
                              DipolarSpheres_Epot_reci_get_structure_modulus
        procedure :: Epot_reci_reInit_structure => DipolarSpheres_Epot_reci_reInit_structure
        procedure, private :: Epot_reci_init_potential => DipolarSpheres_Epot_reci_init_potential
        procedure :: Epot_reci_count_waveVectors => DipolarSpheres_Epot_reci_count_waveVectors
        !>     Reciprocal : delta
        procedure, private :: deltaEpot_reci_move => DipolarSpheres_deltaEpot_reci_move
        procedure, private :: deltaEpot_reci_move_update_structure => &
                              DipolarSpheres_deltaEpot_reci_move_update_structure
        procedure, private :: deltaEpot_reci_rotate => DipolarSpheres_deltaEpot_reci_rotate
        procedure, private :: deltaEpot_reci_rotate_update_structure => &
                              DipolarSpheres_deltaEpot_reci_rotate_update_structure
        procedure, private :: deltaEpot_reci_exchange => DipolarSpheres_deltaEpot_reci_exchange
        !>     Reciprocal : total
        procedure, private :: Epot_reci => DipolarSpheres_Epot_reci
        !>     Self
        procedure, private :: Epot_self_solo => DipolarSpheres_Epot_self_solo
        procedure, private :: Epot_self => DipolarSpheres_Epot_self
        !>     Boundary conditions
        procedure :: Epot_bound_init_totalMoment => DipolarSpheres_Epot_bound_init_totalMoment
        procedure, private :: deltaEpot_bound_rotate => DipolarSpheres_deltaEpot_bound_rotate
        procedure, private :: deltaEpot_bound_rotate_update_totalMoment => &
                              DipolarSpheres_deltaEpot_bound_rotate_update_totalMoment
        procedure, private :: deltaEpot_bound_exchange => DipolarSpheres_deltaEpot_bound_exchange
        procedure, private :: deltaEpot_bound_exchange_update_totalMoment => &
                              DipolarSpheres_deltaEpot_bound_exchange_update_totalMoment
        procedure :: Epot_bound_reInit_totalMoment => DipolarSpheres_Epot_bound_reInit_totalMoment
        procedure, private :: Epot_bound => DipolarSpheres_Epot_bound
        !>     Total
        procedure :: Epot_conf => DipolarSpheres_Epot_conf
        procedure :: test_consist => DipolarSpheres_test_consist
        
    end type DipolarSpheres
    
contains

    subroutine DipolarSpheres_construct(this)
    
        class(DipolarSpheres), intent(out) :: this
        
        this%name = "dipol"        
        write(output_unit, *) this%name, " class construction"
    
        ! Particles
        this%rMin = 1._DP
        this%radius = this%rMin/2._DP
        this%Ncol = dipol_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        allocate(this%orientations(Ndim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = dipol_snap_factor
        
        ! Monte-Carlo
        this%move_delta = dipol_move_delta
        this%move_deltaSave = this%move_delta
        this%move_rejectFix = dipol_move_rejectFix
        this%move_Nadapt = dipol_move_Nadapt
        this%structure_iStep = dipol_structure_iStep
        this%totalMoment_iStep = dipol_totalMoment_iStep
        
        this%rotate_delta = dipol_rotate_delta
        this%rotate_deltaSave = this%rotate_delta
        this%rotate_deltaMax = dipol_rotate_deltaMax
        this%rotate_rejectFix = dipol_rotate_rejectFix
        this%rotate_Nadapt = dipol_rotate_Nadapt
        
        this%Nwidom = dipol_Nwidom
        
        ! Potential
        this%rCut = dipol_rCut
        this%dr = dipol_dr
        this%iMin = int(this%rMin/this%dr)
        this%iCut = int(this%rCut/this%dr)
        this%alpha = dipol_alpha        
        allocate(this%Epot_real_tab(this%iMin:this%iCut, 2))
        call this%Epot_real_init()

        allocate(this%Epot_reci_potential(Ndim, this%Ncol))
        call this%Epot_reci_init_weight()
        
        ! Neighbour Cells
        call this%sameCells%construct(dipol_cell_size, this%rCut) !< same kind
    
    end subroutine DipolarSpheres_construct
    
    subroutine DipolarSpheres_destroy(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
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
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine DipolarSpheres_destroy
    
    !> Report
    
    subroutine DipolarSpheres_print_report(this, report_unit)
    
        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: report_unit    
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        write(report_unit ,*) "    move_Nadapt = ", this%move_Nadapt

        write(report_unit, *) "    alpha = ", this%alpha
        write(report_unit, *) "    rCut = ", this%rCut
        write(report_unit, *) "    dr = ", this%dr
        write(report_unit, *) "    Structure_iStep = ", this%structure_iStep
        write(report_unit, *) "    NwaveVectors = ", this%NwaveVectors
        
        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
        
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
    
    subroutine DipolarSpheres_adapt_rotate_delta(this, reject)
    
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
        
        if (reject == 0._DP) then
            write(error_unit, *) this%name, " :    Warning : rotate_delta adaptation problem."
            this%rotate_delta = this%rotate_deltaSave
            write(error_unit, *) "default rotate_delta :", this%rotate_delta
        end if
        
        if (this%rotate_delta > this%rotate_deltaMax) then
            write(error_unit, *) this%name, " :   Warning : rotate_delta too big."
            this%rotate_delta = this%rotate_deltaMax
            write(error_unit, *) "big rotate_delta :", this%rotate_delta
        end if
        
        write(output_unit, *) this%name, " :    Thermalisation : over (rotation)"
        
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
    
    !> Accessor : move_Nadapt
    
    pure function DipolarSpheres_get_rotate_Nadapt(this) result(get_rotate_Nadapt)
    
        class(DipolarSpheres), intent(in) :: this        
        integer :: get_rotate_Nadapt
        
        get_rotate_Nadapt = this%rotate_Nadapt
        
    end function DipolarSpheres_get_rotate_Nadapt
    
    !> Accessor : structure_iStep
    
    pure function DipolarSpheres_get_structure_iStep(this) result (get_structure_iStep)
    
        class(DipolarSpheres), intent(in) :: this
        integer :: get_structure_iStep
    
        get_structure_iStep = this%structure_iStep
        
    end function DipolarSpheres_get_structure_iStep
    
    !> Accessor : totalMoment_iStep
    
    pure function DipolarSpheres_get_totalMoment_iStep(this) result (get_totalMoment_iStep)
    
        class(DipolarSpheres), intent(in) :: this
        integer :: get_totalMoment_iStep
    
        get_totalMoment_iStep = this%totalMoment_iStep
        
    end function DipolarSpheres_get_totalMoment_iStep

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

    subroutine DipolarSpheres_Epot_real_test_overlap(this, iCol, xCol, iTotalCell, overlap)

        class(DipolarSpheres), intent(in) :: this
        integer, intent(in) :: iCol, iTotalCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap

        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij

        type(Link), pointer :: current => null(), next => null()

        overlap = .false.        

        do iNearCell = 1, NnearCell

            nearCell_index = this%sameCells%nearCells_among_totalCells(iNearCell, iTotalCell)
            current => this%sameCells%beginCells(nearCell_index)%particle%next
            if (.not. associated(current%next)) cycle

            do

                next => current%next

                if (current%iCol /= iCol) then

                    r_ij = dist_PBC(xCol(:), this%positions(:, current%iCol))

                    if (r_ij < this%rMin) then
                        overlap = .true.
                        return
                    end if

                end if

                if (.not. associated(next%next)) exit

                current => next

            end do

        end do

    end subroutine DipolarSpheres_Epot_real_test_overlap
    
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
    
    !> Initialisation of the ``potential'' and the ``structure factor''

    subroutine DipolarSpheres_Epot_reci_init(this)

        class(DipolarSpheres), intent(inout) :: this

        call this%Epot_reci_init_structure()
        call this%Epot_reci_init_potential()

    end subroutine
    
    !> \f[ 
    !>      w(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}{\alpha^2} \sum_d \frac{k_d^2}{L_d}}}
    !>                                {\sum_d \frac{k_d^2}{L_d}}
    !> \f]
    
    subroutine DipolarSpheres_Epot_reci_init_weight(this)
        
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

            if (norm2(waveVector) /= 0._DP) then
            
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

    subroutine DipolarSpheres_Epot_reci_init_structure(this)

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

        this%structureFactor(:, :, :) = cmplx(0._DP, 0._DP, DP)

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
                          
                this%structureFactor(kx, ky, kz) = this%structureFactor(kx, ky, kz) + &
                                                   cmplx(k_dot_mCol, 0._DP, DP) * exp_IkxCol
            
            end do
            
            end do
            
            end do
            
        end do

    end subroutine DipolarSpheres_Epot_reci_init_structure
    
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
                                                  abs(this%structureFactor(kx, ky, kz))

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
    
    !> Potential initialisation :
    !> \f[
    !>      \vec{\phi}(\vec{x}_j) = \sum_{\vec{k}\neq\vec{0}} \vec{k} w(\alpha, \vec{k}) S(\vec{k})
    !>                              e^{-i\vec{k}\cdot\vec{x}_j}
    !> \f]
    
    subroutine DipolarSpheres_Epot_reci_init_potential(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        complex(DP) :: conjg_exp_IkxCol
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_Ikx_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_Ikx_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_Ikx_3
        complex(DP) :: Epot_reci_weight

        complex(DP), dimension(Ndim) :: potential_k
        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: waveVector
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

                Epot_reci_weight = cmplx(this%Epot_reci_weight(kx, ky, kz), 0._DP, DP)

                potential_k(:) = cmplx(waveVector(:), 0._DP, DP) * Epot_reci_weight * &
                                 this%structureFactor(kx, ky, kz)
                
                this%Epot_reci_potential(:, iCol) = this%Epot_reci_potential(:, iCol) + &
                                                    potential_k(:) * conjg_exp_IkxCol

                
            end do
            end do
            end do
            
        end do
        
    end subroutine DipolarSpheres_Epot_reci_init_potential

    ! Count the number of wave vectors

    subroutine DipolarSpheres_Epot_reci_count_waveVectors(this, waveVectors_unit)

        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: waveVectors_unit
        
        real(DP), dimension(Ndim) :: waveVector
        integer :: kx, ky, kz

        this%NwaveVectors = 0

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -Kmax1_sym(ky, kz), Kmax(1)

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

    !> Implementation :
    !> \f[
    !>  \Delta M^2 = 2(\vec{\mu_l}\cdot\vec{k})
    !>              [\cos(\vec{k}\cdot\vec{x}^\prime_l) - \cos(\vec{k}\cdot\vec{x}_l)]
    !>              [\Re(S(\vec{k})) - (\vec{k}\cdot\vec{\mu}_l) \cos(\vec{k}\cdot\vec{x}_l)] -
    !>              [-\sin(\vec{k}\cdot\vec{x}^\prime_l) + \sin(\vec{k}\cdot\vec{x}_l)]
    !>              [\Im(S(\vec{k})) - (\vec{k}\cdot\vec{\mu}_l) \sin(\vec{k}\cdot\vec{x}_l)]
    !> \f]
    !>

    pure function DipolarSpheres_deltaEpot_reci_move(this, xOld, xNew, mCol) &
                  result(deltaEpot_reci_move)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xOld, xNew
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_reci_move
        
        real(DP) :: deltaEpot_k

        real(DP), dimension(Ndim) :: xNewOverL, xOldOverL
        real(DP), dimension(Ndim) :: mColOverL

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
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        complex(DP) :: structure_k
        integer :: kx, ky, kz

        xNewOverL(:) = xNew(:)/Lsize(:)
        xOldOverL(:) = xOld(:)/Lsize(:)
        
        call fourier(xNewOverL, exp_IkxNew_1, exp_IkxNew_2, exp_IkxNew_3)
        call fourier(xOldOverL, exp_IkxOld_1, exp_IkxOld_2, exp_IkxOld_3)

        mColOverL(:) = mCol(:)/Lsize(:)

        deltaEpot_reci_move = 0._DP

        do kz = 0, Kmax(3) ! symmetry : half wave vectors -> double Energy

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2) 

                waveVector(2) = real(ky, DP)
            
                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)

                    k_dot_mCol = dot_product(waveVector, mColOverL)

                    structure_k = this%structureFactor(kx, ky, kz)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    cos_kxNew = real(exp_IkxNew, DP)
                    sin_kxNew = aimag(exp_IkxNew)

                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)
                    cos_kxOld = real(exp_IkxOld, DP)
                    sin_kxOld = aimag(exp_IkxOld)

                    realPart1 = (cos_kxNew - cos_kxOld)
                    realPart1 = realPart1 * (real(structure_k, DP) - k_dot_mCol * cos_kxOld)

                    realPart2 = (-sin_kxNew + sin_kxOld)
                    realPart2 = realPart2 * (aimag(structure_k) - k_dot_mCol * sin_kxOld)
                    
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
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    subroutine DipolarSpheres_deltaEpot_reci_move_update_structure(this, xOld, xNew, mCol)

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

        xNewOverL(:) = xNew(:)/Lsize(:)
        xOldOverL(:) = xOld(:)/Lsize(:)
        
        call fourier(xNewOverL, exp_IkxNew_1, exp_IkxNew_2, exp_IkxNew_3)
        call fourier(xOldOverL, exp_IkxOld_1, exp_IkxOld_2, exp_IkxOld_3)

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
                                                          
                    this%structureFactor(kx, ky, kz) = this%structureFactor(kx, ky, kz) + &
                        cmplx(k_dot_mCol, 0._DP, DP) * (exp_IkxNew - exp_IkxOld)

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_deltaEpot_reci_move_update_structure
    
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
    
    !> Implementation :
    !> \f[
    !>  \Delta M^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2 [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>               \{
    !>                  \cos(\vec{k} \cdot \vec{x}_l)[\Re(S(\vec{k})) -
    !>                      (\vec{k} \cdot \vec{\mu}_l) \cos(\vec{k} \cdot \vec{x}_l)] +
    !>                  \sin(\vec{k} \cdot \vec{x}_l)[\Im(S(\vec{k})) -
    !>                      (\vec{k} \cdot \vec{\mu}_l) \sin(\vec{k} \cdot \vec{x}_l)]
    !>               \}
    !> \f]

    pure function DipolarSpheres_deltaEpot_reci_rotate(this, xCol, mOld, mNew) &
                  result(deltaEpot_reci_rotate)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_reci_rotate

        real(DP) :: deltaEpot_k

        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mNewOverL, mOldOverL

        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol
        real(DP) :: cos_kxCol, sin_kxCol

        real(DP) :: realPart, realPart1, realPart2
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mNew, k_dot_mOld
        complex(DP) :: structure_k
        integer :: kx, ky, kz

        xColOverL(:) = xCol(:)/Lsize(:)
        
        call fourier(xColOverL, exp_IkxCol_1, exp_IkxCol_2, exp_IkxCol_3)

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

                    structure_k = this%structureFactor(kx, ky, kz)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)
                    cos_kxCol = real(exp_IkxCol, DP)
                    sin_kxCol = aimag(exp_IkxCol)

                    realPart1 = cos_kxCol * (real(structure_k, DP)-k_dot_mOld*cos_kxCol)
                    realPart2 = sin_kxCol * (aimag(structure_k)-k_dot_mOld*sin_kxCol)

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
    !>      \Delta S(\vec{k}) = [\vec{k}\cdot(\vec{\mu}_l^\prime - \vec{\mu}_l)] 
    !>                          e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    subroutine DipolarSpheres_deltaEpot_reci_rotate_update_structure(this, xCol, mOld, mNew)

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

        xColOverL(:) = xCol(:)/Lsize(:)
        
        call fourier(xColOverL, exp_IkxCol_1, exp_IkxCol_2, exp_IkxCol_3)

        mNewOverL(:) = mNew(:)/Lsize(:)
        mOldOverL(:) = mOld(:)/Lsize(:)

        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)

                do kx = -Kmax1_sym(ky, kz), Kmax(1)

                    waveVector(1) = real(kx, DP)

                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)

                    k_dot_deltaMcol = dot_product(waveVector, mNewOverL - mOldOverL)

                    this%structureFactor(kx, ky, kz) = this%structureFactor(kx, ky, kz) + &
                        cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_IkxCol

                end do
                
            end do
            
        end do

    end subroutine DipolarSpheres_deltaEpot_reci_rotate_update_structure

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
    
    !> Implementation :
    !> \f[ \Delta U^{N+1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}}
    !>                          (\vec{k} \cdot +\vec{\mu}_{N+1}) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot +\vec{\mu}_{N+1}) +
    !>                              2 [\Re(S(\vec{k})) \cos(\vec{k} \cdot \vec{x}_{N+1}) +
    !>                                 \Im(S(\vec{k})) \sin(\vec{k} \cdot \vec{x}_{N+1})]
    !>                          \}
    !> \f]
    
    !> Subtraction : 
    
    !> Difference of Energy 
    !> \f[ \Delta U_{N-1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}} 
    !>                          (\vec{k} \cdot -\vec{\mu}_N) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot -\vec{\mu}_N) + 
    !>                              2\Re[S_N(\vec{k}) e^{-i \vec{k} \cdot \vec{x}_N}]
    !>                          \}
    !> \f]
    
    !> Implementation :
    !> \f[ \Delta _{N-1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}}
    !>                          (\vec{k} \cdot -\vec{\mu}_N) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot -\vec{\mu}_N) +
    !>                              2 [\Re(S_N(\vec{k})) \cos(\vec{k} \cdot \vec{x}_N) +
    !>                                 \Im(S_N(\vec{k})) \sin(\vec{k} \cdot \vec{x}_N)]
    !>                          \}
    !> \f]
    
    !> Summary : only the sign of \f[\vec{\mu}\f] changes.

    pure function DipolarSpheres_deltaEpot_reci_exchange(this, xCol, mCol) &
                  result(deltaEpot_reci_exchange)

        class(DipolarSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: xCol
        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_reci_exchange
        
        real(DP) :: deltaEpot_k
        
        real(DP), dimension(Ndim) :: xColOverL
        real(DP), dimension(Ndim) :: mColOverL
        
        complex(DP), dimension(-Kmax(1):Kmax(1)) :: exp_IkxCol_1
        complex(DP), dimension(-Kmax(2):Kmax(2)) :: exp_IkxCol_2
        complex(DP), dimension(-Kmax(3):Kmax(3)) :: exp_IkxCol_3
        complex(DP) :: exp_IkxCol
        real(DP) :: cos_kxCol, sin_kxCol
        
        real(DP) :: realPart
        
        real(DP), dimension(Ndim) :: waveVector
        real(DP) :: k_dot_mCol
        complex(DP) :: structure_k
        integer :: kx, ky, kz
        
        xColOverL(:) = xCol(:)/Lsize(:)
        
        call fourier(xColOverL, exp_IkxCol_1, exp_IkxCol_2, exp_IkxCol_3)
        
        mColOverL(:) = mCol(:)/Lsize(:)
        
        deltaEpot_reci_exchange = 0._DP
        
        do kz = 0, Kmax(3)

            waveVector(3) = real(kz, DP)

            do ky = -Kmax2_sym(kz), Kmax(2)

                waveVector(2) = real(ky, DP)
            
                do kx = -Kmax1_sym(ky, kz), Kmax(1)
                
                    waveVector(1) = real(kx, DP)
                    
                    k_dot_mCol = dot_product(waveVector, mColOverL)
                    
                    structure_k = this%structureFactor(kx, ky, kz)
                                                  
                    exp_IkxCol = exp_IkxCol_1(kx) * exp_IkxCol_2(ky) * exp_IkxCol_3(kz)
                    cos_kxCol = real(exp_IkxCol, DP)
                    sin_kxCol = aimag(exp_IkxCol)
                    
                    realPart = real(structure_k, DP) * cos_kxCol
                    realPart = realPart + aimag(structure_k) * sin_kxCol
                    
                    deltaEpot_k = k_dot_mCol * (k_dot_mCol + 2._DP * realPart)
                    deltaEpot_k = deltaEpot_k * this%Epot_reci_weight(kx, ky, kz)
                    deltaEpot_reci_exchange = deltaEpot_reci_exchange + deltaEpot_k
                   
                end do
            
            end do
        
        end do
        
        deltaEpot_reci_exchange = 4._DP*PI/Volume * deltaEpot_reci_exchange

    end function DipolarSpheres_deltaEpot_reci_exchange
    
    !> Total reciprocal energy
    
    pure function DipolarSpheres_Epot_reci(this) result(Epot_reci)
        
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_reci
        
        integer :: jCol
        real(DP), dimension(Ndim) :: mColOverL
        real(DP), dimension(Ndim) :: real_potential

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
    
    ! Boundary conditions : shape-dependent -------------------------------------------------------
    
    !> Total dipole moment :
    !> \f[ \vec{M} = \sum_j \vec{\mu}_j \f]
    !> \f[ \vec{M}_l = \sum_{j \neq l} \vec{\mu}_j \f]
    
    subroutine DipolarSpheres_Epot_bound_init_totalMoment(this)
    
        class(DipolarSpheres), intent(inout) :: this
        
        integer :: iCol
        
        this%totalMoment(:) = 0._DP
        
        do iCol = 1, this%Ncol
        
            this%totalMoment(:) = this%totalMoment(:) + this%orientations(:, iCol)
        
        end do        
    
    end subroutine DipolarSpheres_Epot_bound_init_totalMoment
    
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
                          
        deltaEpot_bound_rotate = 2._DP*PI / (2._DP*out_permittivity+1._DP) / Volume * &
                                 deltaEpot_bound_rotate
    
    end function DipolarSpheres_deltaEpot_bound_rotate
    
    !> Rotation
    
    !> Update the total moment
    !> \f[
    !>      \Delta \vec{M} = \vec{\mu}^\prime_l - \vec{\mu}_l
    !> \f]    
    
    subroutine DipolarSpheres_deltaEpot_bound_rotate_update_totalMoment(this, mOld, mNew)
    
        class(DipolarSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mOld, mNew
        
        this%totalMoment(:) = this%totalMoment(:) + mNew(:) - mOld(:)
    
    end subroutine DipolarSpheres_deltaEpot_bound_rotate_update_totalMoment
    
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
                          
        deltaEpot_bound_exchange = 2._DP*PI / (2._DP*out_permittivity+1._DP) / Volume * &
                                   deltaEpot_bound_exchange
    
    end function DipolarSpheres_deltaEpot_bound_exchange
    
    !> Exchange
    
    !> Update the total moment : add
    !> \f[
    !>      \Delta \vec{M} = +\vec{\mu}_l
    !> \f]
    
    !> Update the total moment : remove
    !> \f[
    !>      \Delta \vec{M} = -\vec{\mu}_l
    !> \f]  
    
    subroutine DipolarSpheres_deltaEpot_bound_exchange_update_totalMoment(this, mCol)
    
        class(DipolarSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mCol
        
        this%totalMoment(:) = this%totalMoment(:) + mCol(:)
    
    end subroutine DipolarSpheres_deltaEpot_bound_exchange_update_totalMoment
    
    !> Reinitialise the total moment factor and print the drift
    
    subroutine DipolarSpheres_Epot_bound_reInit_totalMoment(this, iStep, modulus_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: modulus_unit
        
        real(DP) :: modulus_drifted, modulus_reInit
        
        modulus_drifted = norm2(this%totalMoment(:))
        call this%Epot_bound_init_totalMoment()
        modulus_reInit = norm2(this%totalMoment(:))
        
        write(modulus_unit, *) iStep, abs(modulus_reInit - modulus_drifted)
    
    end subroutine DipolarSpheres_Epot_bound_reInit_totalMoment
    
    !> Total shape dependent term
    !> \f[
    !>      J(\vec{M}, S) = \frac{2\pi}{(2\epsilon_s+1)V} | \vec{M}|^2
    !> \f]
    
    pure function DipolarSpheres_Epot_bound(this) result(Epot_bound)
    
        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_bound
        
        Epot_bound = 2._DP * PI / (2._DP*out_permittivity + 1._DP) / Volume * &
                     dot_product(this%totalMoment, this%totalMoment)
    
    end function DipolarSpheres_Epot_bound

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

