!> \brief Description of the Hard Spheres class

module class_hard_spheres

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use data_box, only: num_dimensions
use json_module, only: json_file
use module_data, only: test_data_found
use module_physics_micro, only: PBC_distance

implicit none

private

    type, public :: Hard_Spheres
        
        private
    
        character(len=:), allocatable :: name

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
        procedure, private :: set_snap => Hard_Spheres_set_snap
        procedure :: destroy => Hard_Spheres_destroy
        
        !> Accessors & Mutators
        procedure :: get_name => Hard_Spheres_get_name
        procedure :: get_num_particles => Hard_Spheres_get_num_particles
        procedure :: get_widom_num_particles => Hard_Spheres_get_widom_num_particles
        procedure :: get_diameter => Hard_Spheres_get_diameter
        procedure :: get_position => Hard_Spheres_get_position
        procedure :: set_position => Hard_Spheres_set_position
        
        procedure :: write_report => Hard_Spheres_write_report
        
        procedure :: write_snap_data => Hard_Spheres_write_snap_data
        procedure :: write_snap_positions => Hard_Spheres_write_snap_positions
        
        procedure :: test_overlap => Hard_Spheres_test_overlap
        
    end type Hard_Spheres
    
    type, extends(Hard_Spheres), public :: Dipolar_Hard_Spheres

        private
        
        ! Particles
        real(DP), dimension(:, :), allocatable, public :: all_orientations
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Dipolar_Hard_Spheres_construct
        procedure, private :: set_particles => Dipolar_Hard_Spheres_set_particles
        procedure :: destroy => Dipolar_Hard_Spheres_destroy
        
        !> Accessor & Mutator
        procedure :: get_orientation => Dipolar_Hard_Spheres_get_orientation
        procedure :: set_orientation => Dipolar_Hard_Spheres_set_orientation
        
        procedure :: write_snap_orientations => Dipolar_Hard_Spheres_write_snap_orientations
        
    end type Dipolar_Hard_Spheres
    
    type, public :: Between_Hard_Spheres
    
        character(len=:), allocatable :: name
        real(DP) :: non_additivity
        real(DP) :: diameter
        
    contains
        
        procedure :: construct => Between_Hard_Spheres_construct
        procedure :: destroy => Between_Hard_Spheres_destroy
        
        procedure :: get_name => Between_Hard_Spheres_get_name
        procedure :: get_diameter => Between_Hard_Spheres_get_diameter
        procedure :: test_overlap => Between_Hard_Spheres_test_overlap
    
    end type Between_Hard_Spheres
    
contains

    subroutine Hard_Spheres_construct(this, json)    
    
        class(Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        character(len=:), allocatable :: this_name
        
        data_name = "Particles.Hard Spheres.name"
        call json%get(data_name, this_name, found)
        call test_data_found(data_name, found)
        this%name = this_name
        if(allocated(this_name)) deallocate(this_name)
        
        write(output_unit, *) this%name, " class construction"
        
        call this%set_particles(json)
        call this%set_snap(json)
        
    end subroutine Hard_Spheres_construct
    
    subroutine Dipolar_Hard_Spheres_construct(this, json)
    
        class(Dipolar_Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        character(len=:), allocatable :: this_name
        
        data_name = "Particles.Dipolar Hard Spheres.name"
        call json%get(data_name, this_name, found)
        call test_data_found(data_name, found)
        this%name = this_name
        if(allocated(this_name)) deallocate(this_name)
        
        write(output_unit, *) this%name, " class construction"
    
        call this%set_particles(json)        
        call this%Hard_Spheres%set_snap(json)
    
    end subroutine Dipolar_Hard_Spheres_construct
    
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
        allocate(this%all_positions(num_dimensions, this%num_particles))
        
        data_name = "Particles.Hard Spheres.number of Widom particles"
        call json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
        
    end subroutine Hard_Spheres_set_particles
    
    subroutine Dipolar_Hard_Spheres_set_particles(this, json)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: json
        
        character(len=4096) :: data_name
        logical :: found
        
        this%diameter = 1._DP ! = u_length
        
        data_name = "Particles.Dipolar Hard Spheres.number of particles"
        call json%get(data_name, this%num_particles, found)
        call test_data_found(data_name, found)
        allocate(this%all_positions(num_dimensions, this%num_particles))
        allocate(this%all_orientations(num_dimensions, this%num_particles))
        
        data_name = "Particles.Dipolar Hard Spheres.number of Widom particles"
        call json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
        
    end subroutine Dipolar_Hard_Spheres_set_particles
    
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
    
    subroutine Between_Hard_Spheres_construct(this, json, type1_diameter, type2_diameter)
    
        class(Between_Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: json
        real(DP), intent(in) :: type1_diameter, type2_diameter
        
        character(len=4096) :: data_name
        logical :: found
        
        character(len=:), allocatable :: this_name
        
        data_name = "Particles.Between Spheres.name"
        call json%get(data_name, this_name, found)
        call test_data_found(data_name, found)
        this%name = this_name
        if(allocated(this_name)) deallocate(this_name)
        
        write(output_unit, *) this%name, " class construction"
        
        data_name = "Particles.Between Spheres.non addivity"
        call json%get(data_name, this%non_additivity, found)
        call test_data_found(data_name, found)        
        this%diameter = (type1_diameter + type2_diameter)/2._DP + this%non_additivity
        
    end subroutine Between_Hard_Spheres_construct
    
    subroutine Hard_Spheres_destroy(this)
    
        class(Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"        
        if (allocated(this%all_positions)) deallocate(this%all_positions)
        if (allocated(this%name)) deallocate(this%name)
    
    end subroutine Hard_Spheres_destroy
    
    subroutine Dipolar_Hard_Spheres_destroy(this)    
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        
        call this%Hard_Spheres%destroy()
        if (allocated(this%all_orientations)) deallocate(this%all_orientations) 
          
    end subroutine Dipolar_Hard_Spheres_destroy
    
    subroutine Between_Hard_Spheres_destroy(this)
    
        class(Between_Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"        
        if (allocated(this%name)) deallocate(this%name)
    
    end subroutine Between_Hard_Spheres_destroy
    
    !> Accessors
    
    pure function Hard_Spheres_get_name(this) result(get_name)
    
        class(Hard_Spheres), intent(in) :: this
        character(len=len(this%name)) :: get_name
        
        get_name = this%name
        
    end function Hard_Spheres_get_name
    
    pure function Between_Hard_Spheres_get_name(this) result(get_name)
    
        class(Between_Hard_Spheres), intent(in) :: this
        character(len=len(this%name)) :: get_name
        
        get_name = this%name
        
    end function Between_Hard_Spheres_get_name

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
        real(DP), dimension(num_dimensions) :: get_position
        
        get_position(:) = this%all_positions(:, i_particle)
        
    end function Hard_Spheres_get_position
    
    subroutine Hard_Spheres_set_position(this, i_particle, position)
    
        class(Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(:), intent(in) :: position
    
        this%all_positions(:, i_particle) = position(:)
        
    end subroutine Hard_Spheres_set_position
    
    pure function Dipolar_Hard_Spheres_get_orientation(this, i_particle) result(get_orientation)
    
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(num_dimensions) :: get_orientation
        
        get_orientation(:) = this%all_orientations(:, i_particle)
        
    end function Dipolar_Hard_Spheres_get_orientation
    
    subroutine Dipolar_Hard_Spheres_set_orientation(this, i_particle, orientation)
    
        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), dimension(:), intent(in) :: orientation
        
        this%all_orientations(:, i_particle) = orientation(:)
        
    end subroutine Dipolar_Hard_Spheres_set_orientation
    
    pure function Between_Hard_Spheres_get_diameter(this) result(get_diameter)
    
        class(Between_Hard_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
        
    end function Between_Hard_Spheres_get_diameter
    
    !> Write a report of the component in a file
    
    subroutine Hard_Spheres_write_report(this, report_unit)
    
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: report_unit
        
        write(report_unit, *) "Data: "
        
        write(report_unit, *) "    snap_factor = ", this%snap_factor
        
    end subroutine Hard_Spheres_write_report
    
    !> Take a snap shot of the configuration
    
    subroutine Hard_Spheres_write_snap_data(this, snap_unit)
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        
        write(snap_unit, *) this%name, this%num_particles, this%snap_factor
        
    end subroutine Hard_Spheres_write_snap_data
      
    subroutine Hard_Spheres_write_snap_positions(this, i_step, snap_unit)
        
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_positions(:, i_particle)
            end do
        end if

    end subroutine Hard_Spheres_write_snap_positions
    
    subroutine Dipolar_Hard_Spheres_write_snap_orientations(this, i_step, snap_unit)
        
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit
    
        integer :: i_particle
        
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num_particles
                write(snap_unit, *) this%all_orientations(:, i_particle)
            end do
        end if

    end subroutine Dipolar_Hard_Spheres_write_snap_orientations
    
    !> Do an overlap test
    
    subroutine Hard_Spheres_test_overlap(this, Box_size)
    
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
    
        integer :: i_particle, j_particle
        real(DP), dimension(num_dimensions) :: position_i, position_j
        real(DP) :: r_ij
    
        do j_particle = 1, this%num_particles
            do i_particle = j_particle+1, this%num_particles
            
                position_i(:) = this%all_positions(:, i_particle)
                position_j(:) = this%all_positions(:, j_particle)
                r_ij = PBC_distance(Box_size, position_i, position_j)
                
                if (r_ij < this%diameter) then
                    write(error_unit, *) this%name, "    Overlap !", i_particle, j_particle
                    write(error_unit, *) "    r_ij = ", r_ij
                    error stop
                end if
                    
            end do
        end do
        
        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Hard_Spheres_test_overlap
    
    subroutine Between_Hard_Spheres_test_overlap(this, Box_size, type1, type2)
    
        class(Between_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: type1, type2
        
        integer :: type1_i_particle, type2_i_particle
        real(DP) :: r_mix
        real(DP), dimension(num_dimensions) :: type1_xCol, type2_xCol
        
        do type1_i_particle = 1, type1%get_num_particles()
            do type2_i_particle = 1, type2%get_num_particles()
                    
                type1_xCol(:) = type1%get_position(type1_i_particle)
                type2_xCol(:) = type2%get_position(type2_i_particle)
                r_mix = PBC_distance(Box_size, type1_xCol, type2_xCol)
                if (r_mix < this%diameter) then
                    write(error_unit, *) this%name, ":    Overlap !", type1_i_particle, type2_i_particle
                    write(error_unit, *) "    r_mix = ", r_mix
                    error stop
                end if

            end do
        end do

        write(output_unit, *) this%name, ":    Overlap test: OK !"
    
    end subroutine Between_Hard_Spheres_test_overlap

end module class_hard_spheres
