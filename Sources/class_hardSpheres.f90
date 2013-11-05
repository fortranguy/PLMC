!> \brief Description of the Hard Spheres class

module class_hardSpheres

use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use data_precisions, only : DP, real_zero
use data_constants, only : PI
use data_box, only : Ndim, Lsize, Volume
use data_particles, only : hard_rMin, hard_Ncol
use data_monteCarlo, only : hard_move_delta, hard_move_rejectFix, hard_Nwidom
use data_neighbourCells, only : NnearCell
use data_distribution, only : hard_snap_factor
use module_physics, only : dist_PBC
use class_neighbourCells

implicit none

private

    type, public :: HardSpheres
    
        ! private
        ! The attributes must be private according to the encapsulation principle.
        ! Nevertheless, it is public for inheritance.
    
        character(len=5) :: name

        ! Particles
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: radius !< radius of a particle
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable, public :: positions !< positions of all particles
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        real(DP), dimension(Ndim) :: move_delta !< displacement
        real(DP), dimension(Ndim) :: move_deltaSave
        real(DP) :: move_rejectFix
        integer :: Nwidom

        ! Potential
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(NeighbourCells), public :: sameCells !< same kind
        type(NeighbourCells) :: mixCells !< other kind
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => HardSpheres_construct
        procedure :: destroy => HardSpheres_destroy
        
        !> Accessors
        procedure :: get_name => HardSpheres_get_name
        procedure :: get_Ncol => HardSpheres_get_Ncol
        procedure :: get_Nwidom => HardSpheres_get_Nwidom
        procedure :: get_rMin => HardSpheres_get_rMin
        procedure :: get_radius => HardSpheres_get_radius
        procedure :: get_rCut => HardSpheres_get_rCut
        procedure :: get_move_delta => HardSpheres_get_move_delta
        !> Specifier        
        !> Adapt the displacement move_delta during thermalisation
        procedure :: adapt_move_delta => HardSpheres_adapt_move_delta
        procedure :: set_move_delta => HardSpheres_set_move_delta        
        
        procedure :: print_density => HardSpheres_print_density
        !> Print a report of the component in a file
        procedure :: print_report => HardSpheres_print_report
        
        !> Take a snap shot of the configuration : positions
        procedure :: snap_positions_data => HardSpheres_snap_positions_data
        procedure :: snap_positions => HardSpheres_snap_positions
        
        !> Do an overlap test
        procedure :: test_overlap => HardSpheres_test_overlap
        
        procedure :: construct_mixCells => HardSpheres_construct_mixCells
        !> Assign all particles to cells
        procedure :: all_cols_to_cells => HardSpheres_all_cols_to_cells
        
        !> Potential energy
        procedure :: Epot_print => HardSpheres_Epot_print
        procedure :: Epot_neighCells => HardSpheres_Epot_neighCells
        procedure :: Epot_conf => HardSpheres_Epot_conf
        procedure :: test_consist => HardSpheres_test_consist
        
    end type HardSpheres
    
contains

    subroutine HardSpheres_construct(this)
    
        class(HardSpheres), intent(out) :: this
        
        real(DP), dimension(Ndim) :: cell_size
        
        this%name = "hardS"
        write(output_unit, *) this%name, " class construction"
    
        ! Particles
        this%rMin = hard_rMin
        this%radius = this%rMin/2._DP
        this%Ncol = hard_Ncol
        allocate(this%positions(Ndim, this%Ncol))
        
        ! Snapshot
        this%snap_factor = hard_snap_factor
        
        ! Monte-Carlo
        this%move_delta = hard_move_delta
        this%move_deltaSave = this%move_delta
        this%move_rejectFix = hard_move_rejectFix
        this%Nwidom = hard_Nwidom
                
        ! Potential
        this%rCut = this%rMin
        
        ! Neighbour Cells
        cell_size(:) = this%rCut
        call this%sameCells%construct(cell_size, this%rCut) !< same kind
    
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%positions)) then
            deallocate(this%positions)
        end if
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Accessor : name

    pure function HardSpheres_get_name(this) result(get_name)
    
        class(HardSpheres), intent(in) :: this        
        character(len=5) :: get_name
        
        get_name = this%name
    
    end function HardSpheres_get_name

    !> Accessor : Ncol

    pure function HardSpheres_get_Ncol(this) result(get_Ncol)
    
        class(HardSpheres), intent(in) :: this        
        integer :: get_Ncol
        
        get_Ncol = this%Ncol
    
    end function HardSpheres_get_Ncol

    !> Accessor : Nwidom

    pure function HardSpheres_get_Nwidom(this) result(get_Nwidom)

        class(HardSpheres), intent(in) :: this
        integer :: get_Nwidom

        get_Nwidom = this%Nwidom

    end function HardSpheres_get_Nwidom
    
    !> Accessor : rMin
    
    pure function HardSpheres_get_rMin(this) result(get_rMin)
    
        class(HardSpheres), intent(in) :: this        
        real(DP) :: get_rMin
        
        get_rMin = this%rMin
    
    end function HardSpheres_get_rMin
    
    !> Accessor : radius
    
    pure function HardSpheres_get_radius(this) result(get_radius)
    
        class(HardSpheres), intent(in) :: this        
        real(DP) :: get_radius
        
        get_radius = this%radius
    
    end function HardSpheres_get_radius
    
    !> Accessor : rCut
    
    pure function HardSpheres_get_rCut(this) result(get_rCut)
    
        class(HardSpheres), intent(in) :: this        
        real(DP) :: get_rCut
        
        get_rCut = this%rCut
    
    end function HardSpheres_get_rCut
    
    pure function HardSpheres_get_move_delta(this) result(get_move_delta)
        
        class(HardSpheres), intent(in) :: this
        real(DP) :: get_move_delta
        
        ! average move_delta of 3 vector components
        get_move_delta = sum(this%move_delta)/size(this%move_delta)
        
    end function HardSpheres_get_move_delta
    
    !> Adaptation of move_delta during the thermalisation
    
    pure subroutine HardSpheres_adapt_move_delta(this, reject)
    
        class(HardSpheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        
        real(DP), parameter :: move_delta_eps = 0.05_DP
        real(DP), parameter :: move_reject_eps = 0.1_DP * move_delta_eps
        real(DP), parameter :: more = 1._DP+move_delta_eps
        real(DP), parameter :: less = 1._DP-move_delta_eps
        
        if (reject < this%move_rejectFix - move_reject_eps) then
        
            this%move_delta(:) = this%move_delta(:) * more
            
            if (norm2(this%move_delta) > norm2(Lsize)) then
                this%move_delta(:) = Lsize(:)
            end if
            
        else if (reject > this%move_rejectFix + move_reject_eps) then
        
            this%move_delta(:) = this%move_delta(:) * less
            
        end if
    
    end subroutine HardSpheres_adapt_move_delta
    
    subroutine HardSpheres_set_move_delta(this, reject, report_unit)
    
        class(HardSpheres), intent(inout) :: this
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit

        if (reject < real_zero) then
            write(error_unit, *) this%name, " :    Warning : move_delta adaptation problem."
            this%move_delta(:) = this%move_deltaSave(:)
            write(error_unit, *) "default move_delta :", this%move_delta(:)
        end if

        if (norm2(this%move_delta) > norm2(Lsize)) then
            write(error_unit, *) this%name, " :   Warning : move_delta too big."
            this%move_delta(:) = Lsize(:)
            write(error_unit, *) "big move_delta :", this%move_delta(:)
        end if

        write(report_unit, *) "Displacement :"
        write(report_unit, *) "    move_delta(:) = ", this%move_delta(:)
        write(report_unit, *) "    rejection relative difference = ", &
                                    abs(reject-this%move_rejectFix)/this%move_rejectFix
    
    end subroutine HardSpheres_set_move_delta
    
    !> Print density and compacity
    
    subroutine HardSpheres_print_density(this, total_Ncol, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: total_Ncol
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity, concentration
        
        density = real(this%Ncol + 1, DP) / Volume ! cheating ? cf. Widom
        compacity = 4._DP/3._DP*PI*this%radius**3 * density
        concentration = real(this%Ncol, DP) / real(total_Ncol, DP)
        
        write(output_unit, *) this%name, " : "
        write(output_unit, *) "    density = ", density
        write(output_unit, *) "    compacity = ", compacity
        write(output_unit, *) "    concentration = ", concentration
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
        write(report_unit, *) "    concentration = ", concentration
    
    end subroutine HardSpheres_print_density
    
    !> Report
    
    subroutine HardSpheres_print_report(this, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data :"
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        
        write(report_unit, *) "    rCut = ", this%rCut
        
        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()        
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
        
        write(report_unit, *) "    snap_factor = ", this%snap_factor
        
    end subroutine HardSpheres_print_report
    
    !> Tag the snapshots
    
    subroutine HardSpheres_snap_positions_data(this, snap_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        
        write(snap_unit, *) this%name, this%Ncol, this%snap_factor
    
    end subroutine HardSpheres_snap_positions_data
    
    !> Configuration state : positions
      
    subroutine HardSpheres_snap_positions(this, iStep, snap_unit)
        
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: iCol
        
        if (modulo(iStep, this%snap_factor) == 0) then
        
            do iCol = 1, this%Ncol
                write(snap_unit, *) this%positions(:, iCol)
            end do
            
        end if            

    end subroutine HardSpheres_snap_positions
    
    !> Overlapt test
    
    subroutine HardSpheres_test_overlap(this)
    
        class(HardSpheres), intent(in) :: this
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, this%Ncol
            do iCol = 1, this%Ncol
                if (iCol /= jCol) then
                    
                    r_ij = dist_PBC(this%positions(:, iCol), this%positions(:, jCol))
                    if (r_ij < this%rMin) then
                        write(error_unit, *) this%name, "    Overlap !", iCol, jCol
                        write(error_unit, *) "    r_ij = ", r_ij
                        stop
                    end if
                    
                end if
            end do
        end do
        
        write(output_unit, *) this%name, " :    Overlap test : OK !"
    
    end subroutine HardSpheres_test_overlap
    
     !> Specifier : mixCells construction
    
    subroutine HardSpheres_construct_mixCells(this, mix_cell_size, mix_rCut)
    
        class(HardSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_rCut
        
        write(output_unit, *) this%name, " : mixCells construction"
        
        call this%mixCells%construct(mix_cell_size, mix_rCut)
    
    end subroutine HardSpheres_construct_mixCells
    
    !> Fill cells with colloids
    
    pure subroutine HardSpheres_all_cols_to_cells(this, other)
    
        class(HardSpheres), intent(inout) :: this
        class(HardSpheres), intent(in) :: other
        
        call this%sameCells%all_cols_to_cells(this%Ncol, this%positions)
        call this%mixCells%all_cols_to_cells(other%Ncol, other%positions)
    
    end subroutine HardSpheres_all_cols_to_cells
    
    !> Print the potential : dummy
    
    subroutine HardSpheres_Epot_print(this, Epot_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        write(Epot_unit, *) this%rCut, 0._DP
    
    end subroutine HardSpheres_Epot_print
    
    subroutine HardSpheres_Epot_neighCells(this, iCol, xCol, iTotalCell, overlap, energ)
        
        class(HardSpheres), intent(in) :: this        
        integer, intent(in) :: iCol, iTotalCell
        real(DP), dimension(:), intent(in) :: xCol
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
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
    
    end subroutine HardSpheres_Epot_neighCells
    
    !> Total potential energy : dummy
    
    pure function HardSpheres_Epot_conf(this) result(Epot_conf)
    
        class(HardSpheres), intent(in) :: this        
        real(DP) :: Epot_conf
    
        Epot_conf = this%Ncol * 0._DP
        
    end function HardSpheres_Epot_conf
    
    !> Consistency test : dummy
    
    subroutine HardSpheres_test_consist(this, Epot, report_unit)
    
        class(HardSpheres), intent(in) :: this
        real(DP), intent(in) :: Epot
        integer, intent(in) :: report_unit
        
        real(DP) :: Epot_conf
        real(DP) :: difference
    
        Epot_conf = this%Epot_conf()
        difference = abs(Epot_conf-Epot)
        
        write(report_unit, *) "Consistency test:"
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        write(report_unit, *) "    absolute difference = ", difference
        
        if (difference > 0._DP) then
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine HardSpheres_test_consist

end module class_hardSpheres
