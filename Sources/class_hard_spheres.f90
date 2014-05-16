!> \brief Description of the Hard Spheres class

module class_hard_spheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP, real_zero
use data_constants, only: PI
use data_box, only: Ndim
use json_module, only: json_file
use data_distribution, only: snap_ratio
use module_types_micro, only: Box_Dimensions, Node, Particle_Index
use module_physics_micro, only: dist_PBC
use module_data, only: test_data_found

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
        integer :: widom_num_particles
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Hard_Spheres_construct
        procedure, private :: set_particles => Hard_Spheres_set_particles
        procedure :: set_snap => Hard_Spheres_set_snap
        procedure :: destroy => Hard_Spheres_destroy
        
        !> Accessors & Mutators
        procedure :: get_name => Hard_Spheres_get_name
        procedure :: get_num_particles => Hard_Spheres_get_num_particles
        procedure :: get_widom_num_particles => Hard_Spheres_get_widom_num_particles
        procedure :: get_diameter => Hard_Spheres_get_diameter
        procedure :: get_position => Hard_Spheres_get_position
        procedure :: set_position => Hard_Spheres_set_position    
        
        procedure :: write_density => Hard_Spheres_write_density
        procedure :: write_report => Hard_Spheres_write_report
        
        procedure :: write_snap_data => Hard_Spheres_write_snap_data
        procedure :: write_snap_positions => Hard_Spheres_write_snap_positions
        
        procedure :: test_overlap => Hard_Spheres_test_overlap
        
    end type Hard_Spheres
    
contains

    subroutine Hard_Spheres_construct(this, json)    
        class(Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        this%name = "hardS"
        write(output_unit, *) this%name, " class construction"
        
        call this%set_particles(json)
        call this%set_snap(json)
        
    end subroutine Hard_Spheres_construct
    
    subroutine Hard_Spheres_set_particles(this, json)
        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        data_name = "Particles.Hard Spheres.diameter"
        call json%get(data_name, this%diameter, found)
        call test_data_found(data_name, found)
        
        data_name = "Particles.Hard Spheres.number of particles"
        call json%get(data_name, this%num_particles, found)
        call test_data_found(data_name, found)        
        allocate(this%all_positions(Ndim, this%num_particles))
        
        data_name = "Particles.Hard Spheres.number of Widom particles"
        call json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
        
    end subroutine Hard_Spheres_set_particles
    
    subroutine Hard_Spheres_set_snap(this, json)
    
        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        integer :: snap_ratio
        
        data_name = "Distribution.number of particles"
        call json%get(data_name, snap_ratio, found)
        call test_data_found(data_name, found)
        this%snap_factor = this%num_particles/snap_ratio
        if (this%snap_factor == 0) this%snap_factor = 1
    
    end subroutine Hard_Spheres_set_snap
    
    subroutine Hard_Spheres_destroy(this)
    
        class(Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        
        if (allocated(this%all_positions)) deallocate(this%all_positions)
    
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

    pure function Hard_Spheres_get_widom_num_particles(this) result(get_widom_num_particles)
        class(Hard_Spheres), intent(in) :: this
        integer :: get_widom_num_particles
        
        get_widom_num_particles = this%widom_num_particles
    end function Hard_Spheres_get_widom_num_particles
    
    pure function Hard_Spheres_get_diameter(this) result(get_diameter)
        class(Hard_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
    end function Hard_Spheres_get_diameter
    
    pure function Hard_Spheres_get_position(this, i_particle) result(get_position)
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(Ndim) :: get_position
        
        get_position(:) = this%all_positions(:, i_particle)
    end function Hard_Spheres_get_position
    
    subroutine Hard_Spheres_set_position(this, i_particle, position)
        class(Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(:), intent(in) :: position
    
        this%all_positions(:, i_particle) = position(:)
    end subroutine Hard_Spheres_set_position
    
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
        write(report_unit ,*) "    widom_num_particles = ", this%widom_num_particles
        
        write(report_unit, *) "    snap_factor = ", this%snap_factor
        
    end subroutine Hard_Spheres_write_report
    
    !> Take a snap shot of the configuration: positions
    
    !> Tag the snapshots
    
    subroutine Hard_Spheres_write_snap_data(this, snap_unit)
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        
        write(snap_unit, *) this%name, this%num_particles, this%snap_factor
        
    end subroutine Hard_Spheres_write_snap_data
    
    !> Configuration state: positions
      
    subroutine Hard_Spheres_write_snap_positions(this, iStep, snap_unit)
        
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: iStep
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(iStep, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_positions(:, i_particle)
            end do
        end if

    end subroutine Hard_Spheres_write_snap_positions
    
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
                
                if (r_ij < this%diameter) then
                    write(error_unit, *) this%name, "    Overlap !", i_particle, j_particle
                    write(error_unit, *) "    r_ij = ", r_ij
                    error stop
                end if
                    
            end do
        end do
        
        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Hard_Spheres_test_overlap

end module class_hard_spheres
