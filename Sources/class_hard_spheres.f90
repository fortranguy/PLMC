!> \brief Description of the Hard Spheres class

module class_hard_spheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_constants, only: PI
use data_box, only: Ndim
use data_particles, only: hard_diameter, hard_num_particles
use data_monte_carlo, only: hard_move_delta, hard_move_rejectFix, hard_Nwidom
use data_potential, only: hard_rMin_factor
use data_neighbour_cells, only: NnearCell
use data_distribution, only: snap_ratio
use module_types, only: Box_Dimensions, Node, Particle_Index
use module_physics_micro, only: dist_PBC
use class_neighbour_cells
use class_small_move

implicit none

private

    type, public :: Hard_Spheres
    
        ! private
        ! The attributes must be private according to the encapsulation principle.
        ! Nevertheless, it is public for inheritance.
    
        character(len=5) :: name

        ! Particles
        real(DP) :: diameter
        integer ::  num_particles
        real(DP), dimension(:, :), allocatable, public :: all_positions
        
        ! Snashot
        integer :: snap_factor

        ! Monte-Carlo
        type(Small_Move) :: move
        integer :: Nwidom

        ! Potential
        real(DP) :: rMin !< minimum distance between two particles
        real(DP) :: rCut !< short-range cut
        
        ! Neighbours (cell/grid scheme)
        type(Neighbour_Cells), public :: sameCells !< same kind
        type(Neighbour_Cells) :: mixCells !< other kind
        
    contains

        !> Construction and destruction of the class
        procedure, private :: set_particles => Hard_Spheres_set_particles
        procedure, private :: set_changes => Hard_Spheres_set_changes
        procedure :: construct => Hard_Spheres_construct
        procedure :: destroy => Hard_Spheres_destroy
        
        !> Accessors & Mutators
        procedure :: get_name => Hard_Spheres_get_name
        procedure :: get_num_particles => Hard_Spheres_get_num_particles
        procedure :: get_Nwidom => Hard_Spheres_get_Nwidom
        procedure :: get_diameter => Hard_Spheres_get_diameter
        procedure :: get_move_delta => Hard_Spheres_get_move_delta
        procedure :: get_move_delta_scalar => Hard_Spheres_get_move_delta_scalar
        procedure :: adapt_move_delta => Hard_Spheres_adapt_move_delta
        procedure :: set_move_delta => Hard_Spheres_set_move_delta        
        
        procedure :: write_density => Hard_Spheres_write_density
        procedure :: write_report => Hard_Spheres_write_report
        
        procedure :: snap_data => Hard_Spheres_snap_data
        procedure :: snap_positions => Hard_Spheres_snap_positions
        
        procedure :: test_overlap => Hard_Spheres_test_overlap
        
        !> Neighbour cells
        procedure :: construct_cells => Hard_Spheres_construct_cells
        
        !> Potential energy
        procedure :: set_Epot => Hard_Spheres_set_Epot
        procedure :: write_Epot => Hard_Spheres_write_Epot
        procedure :: Epot_neighCells => Hard_Spheres_Epot_neighCells
        procedure :: Epot_conf => Hard_Spheres_Epot_conf
        
    end type Hard_Spheres
    
contains

    pure subroutine Hard_Spheres_set_particles(this)
        class(Hard_Spheres), intent(inout) :: this
        this%diameter = hard_diameter
        this%num_particles = hard_num_particles
        allocate(this%all_positions(Ndim, this%num_particles))
    end subroutine Hard_Spheres_set_particles
    
    pure subroutine Hard_Spheres_set_changes(this)
        class(Hard_Spheres), intent(inout) :: this
        call this%move%init(hard_move_delta, hard_move_rejectFix)
    end subroutine Hard_Spheres_set_changes

    subroutine Hard_Spheres_construct(this)
    
        class(Hard_Spheres), intent(out) :: this
        
        this%name = "hardS"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_particles()
        call this%set_changes()
        this%Nwidom = hard_Nwidom
        this%snap_factor = this%num_particles/snap_ratio
        if (this%snap_factor == 0) this%snap_factor = 1
        
    end subroutine Hard_Spheres_construct
    
    subroutine Hard_Spheres_destroy(this)
    
        class(Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%all_positions)) deallocate(this%all_positions)
        
        call this%sameCells%destroy()
        call this%mixCells%destroy()
    
    end subroutine Hard_Spheres_destroy
    
    !> Accessors

    pure function Hard_Spheres_get_name(this) result(get_name)
        class(Hard_Spheres), intent(in) :: this
        character(len=5) :: get_name
        
        get_name = this%name
    end function Hard_Spheres_get_name

    pure function Hard_Spheres_get_num_particles(this) result(get_num_particles)
        class(Hard_Spheres), intent(in) :: this
        integer :: get_num_particles
        
        get_num_particles = this%num_particles
    end function Hard_Spheres_get_num_particles

    pure function Hard_Spheres_get_Nwidom(this) result(get_Nwidom)
        class(Hard_Spheres), intent(in) :: this
        integer :: get_Nwidom
        
        get_Nwidom = this%Nwidom
    end function Hard_Spheres_get_Nwidom
    
    pure function Hard_Spheres_get_diameter(this) result(get_diameter)
        class(Hard_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Hard_Spheres_get_diameter
    
    pure function Hard_Spheres_get_move_delta(this) result(get_move_delta)
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(Ndim) :: get_move_delta
        
        get_move_delta = this%move%delta
    end function Hard_Spheres_get_move_delta
    
    pure function Hard_Spheres_get_move_delta_scalar(this) result(get_move_delta_scalar)
        class(Hard_Spheres), intent(in) :: this
        real(DP) :: get_move_delta_scalar
        
        get_move_delta_scalar = sum(this%move%delta)/size(this%move%delta)
    end function Hard_Spheres_get_move_delta_scalar
    
    subroutine Hard_Spheres_adapt_move_delta(this, Box_size, reject)
        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject        
    
        call this%move%adapt_delta(Box_size, reject)
    end subroutine Hard_Spheres_adapt_move_delta
    
    subroutine Hard_Spheres_set_move_delta(this, Box_size, reject, report_unit)
        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        real(DP), intent(in) :: reject
        integer, intent(in) :: report_unit
        
        call this%move%set_delta(this%name, Box_size, reject, report_unit)
    end subroutine Hard_Spheres_set_move_delta
    
    !> Write density and compacity
    
    subroutine Hard_Spheres_write_density(this, Box_size, total_num_particles, report_unit)
    
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size ! warning: average ?
        integer, intent(in) :: total_num_particles
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity, concentration
        
        density = real(this%num_particles + 1, DP) / product(Box_size) ! cheating ? cf. Widom
        compacity = 4._DP/3._DP*PI*(this%diameter/2._DP)**3 * density
        concentration = real(this%num_particles, DP) / real(total_num_particles, DP)
        
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
        write(report_unit, *) "    concentration = ", concentration
    
    end subroutine Hard_Spheres_write_density
    
    !> Write a report of the component in a file
    
    subroutine Hard_Spheres_write_report(this, report_unit)
    
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit ,*) "    num_particles = ", this%num_particles
        write(report_unit ,*) "    Nwidom = ", this%Nwidom
        
        write(report_unit, *) "    rCut = ", this%rCut
        
        write(report_unit, *) "    this_NtotalCell_dim(:) = ", this%sameCells%get_NtotalCell_dim()
        write(report_unit, *) "    this_cell_size(:) = ", this%sameCells%get_cell_size()
        write(report_unit, *) "    mix_NtotalCell_dim(:) = ", this%mixCells%get_NtotalCell_dim()
        write(report_unit, *) "    mix_cell_size(:) = ", this%mixCells%get_cell_size()
        
        write(report_unit, *) "    snap_factor = ", this%snap_factor
        
    end subroutine Hard_Spheres_write_report
    
    !> Take a snap shot of the configuration: positions
    
    !> Tag the snapshots
    
    subroutine Hard_Spheres_snap_data(this, snap_unit)
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        write(snap_unit, *) this%name, this%num_particles, this%snap_factor
    end subroutine Hard_Spheres_snap_data
    
    !> Configuration state: positions
      
    subroutine Hard_Spheres_snap_positions(this, iStep, snap_unit)
        
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(iStep, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_positions(:, i_particle)
            end do
        end if

    end subroutine Hard_Spheres_snap_positions
    
    !> Do an overlap test
    
    subroutine Hard_Spheres_test_overlap(this, Box_size)
    
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
    
        integer :: i_particle, j_particle
        real(DP), dimension(Ndim) :: position_i, position_j
        real(DP) :: r_ij
    
        do j_particle = 1, this%num_particles
            do i_particle = j_particle+1, this%num_particles
            
                position_i(:) = this%all_positions(:, i_particle)
                position_j(:) = this%all_positions(:, j_particle)
                r_ij = dist_PBC(Box_size, position_i, position_j)
                
                if (r_ij < this%rMin) then
                    write(error_unit, *) this%name, "    Overlap !", i_particle, j_particle
                    write(error_unit, *) "    r_ij = ", r_ij
                    error stop
                end if
                    
            end do
        end do
        
        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Hard_Spheres_test_overlap
    
    !> Neighbour Cells
    
    subroutine Hard_Spheres_construct_cells(this, Box_size, other, mix_cell_size, mix_rCut)
    
        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: other
        real(DP), dimension(:), intent(in) :: mix_cell_size
        real(DP), intent(in) :: mix_rCut
        
        real(DP), dimension(Ndim) :: same_cell_size
        
        same_cell_size(:) = this%rCut
        call this%sameCells%construct(Box_size, same_cell_size, this%rCut) !< same kind
        call this%sameCells%all_cols_to_cells(this%num_particles, this%all_positions)
        
        call this%mixCells%construct(Box_size, mix_cell_size, mix_rCut)
        call this%mixCells%all_cols_to_cells(other%num_particles, other%all_positions)
    
    end subroutine Hard_Spheres_construct_cells
    
    ! Potential
    
    subroutine Hard_Spheres_set_Epot(this, Box)
    
        class(Hard_Spheres), intent(inout) :: this
        type(Box_Dimensions), intent(in) :: Box
        
        real(DP) :: volume_dummy
        volume_dummy = product(Box%size)
        
        this%rMin = hard_rMin_factor * this%diameter
        this%rCut = this%rMin
        
    end subroutine Hard_Spheres_set_Epot
    
    !> Write the potential: dummy
    
    subroutine Hard_Spheres_write_Epot(this, Epot_unit)
    
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: Epot_unit

        write(Epot_unit, *) this%rCut, 0._DP
    
    end subroutine Hard_Spheres_write_Epot
    
    subroutine Hard_Spheres_Epot_neighCells(this, Box_size, particle, overlap, energ)
        
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        type(Particle_Index), intent(in) :: particle
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
            
                if (current%number /= particle%number) then
                    r_ij = dist_PBC(Box_size, particle%xCol(:), this%all_positions(:, current%number))
                    if (r_ij < this%rMin) then
                        overlap = .true.
                        return
                    end if
                end if
                
                if (.not. associated(next%next)) exit
                
                current => next
            
            end do
            
        end do
    
    end subroutine Hard_Spheres_Epot_neighCells
    
    !> Total potential energy: dummy
    
    pure function Hard_Spheres_Epot_conf(this, Box) result(Epot_conf)
    
        class(Hard_Spheres), intent(in) :: this
        type(Box_Dimensions), intent(in) :: Box
        real(DP) :: Epot_conf
        
        real(DP) :: volume_dummy
        volume_dummy = product(Box%size)
    
        Epot_conf = this%num_particles * 0._DP
        
    end function Hard_Spheres_Epot_conf

end module class_hard_spheres
