!> \brief Description of the Hard Spheres class

module class_hardSpheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_constants, only: PI
use data_box, only: Ndim
use data_particles, only: hard_sigma, hard_Ncol
use data_monteCarlo, only: hard_move_delta, hard_move_rejectFix, hard_Nwidom
use data_potential, only: hard_rMin_factor
use data_neighbourCells, only: NnearCell
use data_distribution, only: snap_ratio
use module_types, only: Box_dimensions, Node, particle_index
use module_physics_micro, only: dist_PBC
use class_neighbourCells
use class_small_move

implicit none

private

    type, public :: HardSpheres
    
        ! private
        ! The attributes must be private according to the encapsulation principle.
        ! Nevertheless, it is public for inheritance.
    
        character(len=5) :: name

        ! Particles
        real(DP) :: sigma !< diameter
        integer ::  Ncol !< number of a component particles
        real(DP), dimension(:, :), allocatable, public :: positions !< positions of all particles
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        type(Small_move) :: move
        integer :: Nwidom

        ! Potential
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(NeighbourCells), public :: sameCells !< same kind
        type(NeighbourCells) :: mixCells !< other kind
        
    contains

        !> Construction and destruction of the class
        procedure, private :: set_particles => HardSpheres_set_particles
        procedure, private :: set_changes => HardSpheres_set_changes
        procedure :: construct => HardSpheres_construct
        procedure :: destroy => HardSpheres_destroy
        
        !> Accessors & Mutators
        procedure :: get_name => HardSpheres_get_name
        procedure :: get_Ncol => HardSpheres_get_Ncol
        procedure :: get_Nwidom => HardSpheres_get_Nwidom
        procedure :: get_sigma => HardSpheres_get_sigma
        procedure :: get_move_delta => HardSpheres_get_move_delta
        procedure :: adapt_move_delta => HardSpheres_adapt_move_delta
        procedure :: set_move_delta => HardSpheres_set_move_delta        
        
        procedure :: write_density => HardSpheres_write_density
        procedure :: write_report => HardSpheres_write_report
        
        procedure :: snap_data => HardSpheres_snap_data
        procedure :: snap_positions => HardSpheres_snap_positions
        
        procedure :: test_overlap => HardSpheres_test_overlap
        
        !> Neighbour cells
        procedure :: construct_cells => HardSpheres_construct_cells
        
        !> Potential energy
        procedure :: set_Epot => HardSpheres_set_Epot
        procedure :: write_Epot => HardSpheres_write_Epot
        procedure :: Epot_neighCells => HardSpheres_Epot_neighCells
        procedure :: Epot_conf => HardSpheres_Epot_conf
        
    end type HardSpheres
    
contains

    pure subroutine HardSpheres_set_particles(this)
        class(HardSpheres), intent(inout) :: this
        this%sigma = hard_sigma
        this%Ncol = hard_Ncol
        allocate(this%positions(Ndim, this%Ncol))
    end subroutine HardSpheres_set_particles
    
    pure subroutine HardSpheres_set_changes(this)
        class(HardSpheres), intent(inout) :: this
        call this%move%init(hard_move_delta, hard_move_rejectFix)
    end subroutine HardSpheres_set_changes

    subroutine HardSpheres_construct(this)
    
        class(HardSpheres), intent(out) :: this
        
        this%name = "hardS"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_particles()
        call this%set_changes()
        this%Nwidom = hard_Nwidom
        this%snap_factor = this%Ncol/snap_ratio
        if (this%snap_factor == 0) this%snap_factor = 1
        
    end subroutine HardSpheres_construct
    
    subroutine HardSpheres_destroy(this)
    
        class(HardSpheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%positions)) deallocate(this%positions)
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine HardSpheres_destroy
    
    !> Accessors

    pure function HardSpheres_get_name(this) result(get_name)
        class(HardSpheres), intent(in) :: this
        character(len=5) :: get_name
        
        get_name = this%name
    end function HardSpheres_get_name

    pure function HardSpheres_get_Ncol(this) result(get_Ncol)
        class(HardSpheres), intent(in) :: this
        integer :: get_Ncol
        
        get_Ncol = this%Ncol
    end function HardSpheres_get_Ncol

    pure function HardSpheres_get_Nwidom(this) result(get_Nwidom)
        class(HardSpheres), intent(in) :: this
        integer :: get_Nwidom
        
        get_Nwidom = this%Nwidom
    end function HardSpheres_get_Nwidom
    
    pure function HardSpheres_get_sigma(this) result(get_sigma)
        class(HardSpheres), intent(in) :: this
        real(DP) :: get_sigma
        
        get_sigma = this%sigma
    end function HardSpheres_get_sigma
    
    pure function HardSpheres_get_move_delta(this) result(get_move_delta)
        class(HardSpheres), intent(in) :: this
        real(DP), dimension(Ndim) :: get_move_delta
        
        get_move_delta = this%move%delta
    end function HardSpheres_get_move_delta  
    
    subroutine HardSpheres_adapt_move_delta(this, Box_size, reject)
        class(HardSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject        
    
        call this%move%adapt_delta(Box_size, reject)
    end subroutine HardSpheres_adapt_move_delta
    
    subroutine HardSpheres_set_move_delta(this, Box_size, reject, report_unit)
        class(HardSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        call this%move%set_delta(this%name, Box_size, reject, report_unit)
    end subroutine HardSpheres_set_move_delta
    
    !> Write density and compacity
    
    subroutine HardSpheres_write_density(this, Box_size, total_Ncol, report_unit)
    
        class(HardSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average ?
        integer, intent(in) :: total_Ncol
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity, concentration
        
        density = real(this%Ncol + 1, DP) / product(Box_size) ! cheating ? cf. Widom
        compacity = 4._DP/3._DP*PI*(this%sigma/2._DP)**3 * density
        concentration = real(this%Ncol, DP) / real(total_Ncol, DP)
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
        write(report_unit, *) "    concentration = ", concentration
    
    end subroutine HardSpheres_write_density
    
    !> Write a report of the component in a file
    
    subroutine HardSpheres_write_report(this, report_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit ,*) "    Ncol = ", this%Ncol
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        
        write(report_unit, *) "    rCut = ", this%rCut
        
        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
        
        write(report_unit, *) "    snap_factor = ", this%snap_factor
        
    end subroutine HardSpheres_write_report
    
    !> Take a snap shot of the configuration: positions
    
    !> Tag the snapshots
    
    subroutine HardSpheres_snap_data(this, snap_unit)
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        write(snap_unit, *) this%name, this%Ncol, this%snap_factor
    end subroutine HardSpheres_snap_data
    
    !> Configuration state: positions
      
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
    
    !> Do an overlap test
    
    subroutine HardSpheres_test_overlap(this, Box_size)
    
        class(HardSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
    
        integer :: jCol, iCol
        real(DP) :: r_ij
    
        do jCol = 1, this%Ncol
            do iCol = jCol+1, this%Ncol
                    
                r_ij = dist_PBC(Box_size, this%positions(:, iCol), this%positions(:, jCol))
                if (r_ij < this%rMin) then
                    write(error_unit, *) this%name, "    Overlap !", iCol, jCol
                    write(error_unit, *) "    r_ij = ", r_ij
                    error stop
                end if
                    
            end do
        end do
        
        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine HardSpheres_test_overlap
    
    !> Neighbour Cells
    
    subroutine HardSpheres_construct_cells(this, Box_size, other, mix_cell_size, mix_rCut)
    
        class(HardSpheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(HardSpheres), intent(in) :: other
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_rCut
        
        real(DP), dimension(Ndim) :: same_cell_size
        
        same_cell_size(:) = this%rCut
        call this%sameCells%construct(Box_size, same_cell_size, this%rCut) !< same kind
        call this%sameCells%all_cols_to_cells(this%Ncol, this%positions)
        
        call this%mixCells%construct(Box_size, mix_cell_size, mix_rCut)
        call this%mixCells%all_cols_to_cells(other%Ncol, other%positions)
    
    end subroutine HardSpheres_construct_cells
    
    ! Potential
    
    subroutine HardSpheres_set_Epot(this, Box)
    
        class(HardSpheres), intent(inout) :: this
        type(Box_dimensions), intent(in) :: Box
        
        real(DP) :: volume_dummy
        volume_dummy = product(Box%size)
        
        this%rMin = hard_rMin_factor * this%sigma
        this%rCut = this%rMin
        
    end subroutine HardSpheres_set_Epot
    
    !> Write the potential: dummy
    
    subroutine HardSpheres_write_Epot(this, Epot_unit)
    
        class(HardSpheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        write(Epot_unit, *) this%rCut, 0._DP
    
    end subroutine HardSpheres_write_Epot
    
    subroutine HardSpheres_Epot_neighCells(this, Box_size, particle, overlap, energ)
        
        class(HardSpheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(particle_index), intent(in) :: particle
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energ
    
        integer :: iNearCell,  nearCell_index
        real(DP) :: r_ij
    
        type(Node), pointer :: current => null(), next => null()
        
        overlap = .false.
        energ = 0._DP
    
        do iNearCell = 1, NnearCell
        
            nearCell_index = this%sameCells%near_among_total(iNearCell, particle%same_iCell)
            current => this%sameCells%beginCells(nearCell_index)%particle%next
            if (.not. associated(current%next)) cycle
            
            do
            
                next => current%next
            
                if (current%iCol /= particle%iCol) then
                    r_ij = dist_PBC(Box_size, particle%xCol(:), this%positions(:, current%iCol))
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
    
    !> Total potential energy: dummy
    
    pure function HardSpheres_Epot_conf(this, Box) result(Epot_conf)
    
        class(HardSpheres), intent(in) :: this
        type(Box_dimensions), intent(in) :: Box
        real(DP) :: Epot_conf
        
        real(DP) :: volume_dummy
        volume_dummy = product(Box%size)
    
        Epot_conf = this%Ncol * 0._DP
        
    end function HardSpheres_Epot_conf

end module class_hardSpheres
