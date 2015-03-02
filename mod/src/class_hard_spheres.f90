!> \brief Description of the Hard Spheres class

module class_hard_spheres

, output_unit, error_unit
use data_constants, only: PI
use data_box, only: num_dimensions
, json_value, json_value_create, to_object, json_value_add
use data_write, only: simple_precision_format, double_precision_format
use module_data, only: test_data_found, test_empty_string
use module_types_micro, only: Box_Parameters
use module_geometry, only: geometry
use module_physics_micro, only: PBC_distance

implicit none

private

    type, public :: Hard_Spheres
        
        private
    
        

        ! Particles
        
        real(DP) :: volume
        
        ! Snashot
        integer :: snap_factor

        ! Post-Processing
        integer :: widom_num_particles
        
    contains

        !> Construction and destruction of the class
        procedure :: construct => Hard_Spheres_construct
        procedure, private :: set_particles => Hard_Spheres_set_particles
        procedure, private :: set_volume => Hard_Spheres_set_volume
        procedure, private :: set_snap => Hard_Spheres_set_snap
        procedure :: destroy => Hard_Spheres_destroy
        
        !> Accessors & Mutators
        procedure :: get_widom_num_particles => Hard_Spheres_get_widom_num_particles
        procedure :: set_data => Hard_Spheres_set_data
        procedure :: set_positions => Hard_Spheres_set_positions
        procedure :: set_test_num_particles => Hard_Spheres_set_test_num_particles
        procedure, private :: set_widom_num_particles => Hard_Spheres_set_widom_num_particles
        
        procedure :: write_report => Hard_Spheres_write_report
        procedure :: write_data => Hard_Spheres_write_data
        procedure :: write_positions => Hard_Spheres_write_positions
        
        procedure :: test_overlap => Hard_Spheres_test_overlap
        
    end type Hard_Spheres
    
    type, extends(Hard_Spheres), public :: Dipolar_Hard_Spheres

        private
        
        ! Particles
        real(DP), dimension(:, :), allocatable, public :: moments
        
        ! Post-Processing
        integer :: field_num_particles
        
    contains
        
        !> Accessor & Mutator
        procedure :: get_field_num_particles => Dipolar_Hard_Spheres_get_field_num_particles
        procedure :: set_moments => Dipolar_Hard_Spheres_set_moments
        procedure, private :: set_field_num_particles => Dipolar_Hard_Spheres_set_field_num_particles
        
        procedure :: write_moments => Dipolar_Hard_Spheres_write_moments
        
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

    subroutine Hard_Spheres_construct(this, Box, input_data, object_field)
    
        class(Hard_Spheres), intent(out) :: this
        type(Box_Parameters), intent(in) :: Box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: object_field
        
        
        write(output_unit, *) this%name, " class construction"
        
        call this%set_particles(input_data, object_field)
        call this%set_volume(Box)
        call this%set_snap(input_data)
        
    end subroutine Hard_Spheres_construct
    
    subroutine Hard_Spheres_set_particles(this, input_data, object_field)
    
        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: object_field
        
        

        select type (this)
            type is (Hard_Spheres)
                
            type is (Dipolar_Hard_Spheres)
                this%diameter = 1._DP ! = u_length
                allocate(this%moments(num_dimensions, this%num))
        end select
        
    end subroutine Hard_Spheres_set_particles

    subroutine Hard_Spheres_set_volume(this, Box)

        class(Hard_Spheres), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box

        if (geometry%bulk) then
            this%volume = product(Box%size)
        else if (geometry%slab) then
            this%volume = product(Box%size(1:2)) * (Box%height - this%diameter)
        end if

    end subroutine Hard_Spheres_set_volume
    
    subroutine Hard_Spheres_set_snap(this, input_data)
    
        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: input_data
        
        character(len=4096) :: data_field
        logical :: found
        integer :: snap_ratio
        
        data_field = "Distribution.number of particles"
        call input_data%get(data_field, snap_ratio, found)
        call test_data_found(data_field, found)
        this%snap_factor = this%num/snap_ratio
        if (this%snap_factor == 0) this%snap_factor = 1
    
    end subroutine Hard_Spheres_set_snap
    
    subroutine Between_Hard_Spheres_construct(this, input_data, type1_diameter, type2_diameter)
    
        class(Between_Hard_Spheres), intent(out) :: this
        type(json_file), intent(inout) :: input_data
        real(DP), intent(in) :: type1_diameter, type2_diameter
        
        character(len=4096) :: data_field
        logical :: found
        character(len=:), allocatable :: this_name
        
        data_field = "Particles.Between Spheres.name"
        call input_data%get(data_field, this_name, found)
        call test_data_found(data_field, found)
        call test_empty_string(data_field, this_name)
        this%name = this_name
        if (allocated(this_name)) deallocate(this_name)
        write(output_unit, *) this%name, " class construction"
        
        data_field = "Particles.Between Spheres.non addivity"
        call input_data%get(data_field, this%non_additivity, found)
        call test_data_found(data_field, found)
        this%diameter = (type1_diameter + type2_diameter)/2._DP + this%non_additivity
        
    end subroutine Between_Hard_Spheres_construct
    
    subroutine Hard_Spheres_destroy(this)
    
        class(Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        if (allocated(this%positions)) deallocate(this%positions)
        if (allocated(this%name)) deallocate(this%name)

        select type (this)
            type is (Dipolar_Hard_Spheres)
                if (allocated(this%moments)) deallocate(this%moments)
        end select
    
    end subroutine Hard_Spheres_destroy
    
    subroutine Between_Hard_Spheres_destroy(this)
    
        class(Between_Hard_Spheres), intent(inout) :: this
        
        write(output_unit, *) this%name, " class destruction"
        if (allocated(this%name)) deallocate(this%name)
    
    end subroutine Between_Hard_Spheres_destroy
    
    !> Accessors
    
    pure function Between_Hard_Spheres_get_name(this) result(get_name)
    
        class(Between_Hard_Spheres), intent(in) :: this
        character(len=len(this%name)) :: get_name
        
        get_name = this%name
        
    end function Between_Hard_Spheres_get_name

    pure function Hard_Spheres_get_widom_num_particles(this) result(get_widom_num_particles)
    
        class(Hard_Spheres), intent(in) :: this
        integer :: get_widom_num_particles
        
        get_widom_num_particles = this%widom_num_particles
        
    end function Hard_Spheres_get_widom_num_particles

    pure function Dipolar_Hard_Spheres_get_field_num_particles(this) result(get_field_num_particles)

        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer :: get_field_num_particles

        get_field_num_particles = this%field_num_particles

    end function Dipolar_Hard_Spheres_get_field_num_particles
    
    pure function Between_Hard_Spheres_get_diameter(this) result(get_diameter)
    
        class(Between_Hard_Spheres), intent(in) :: this
        real(DP) :: get_diameter
        
        get_diameter = this%diameter
        
    end function Between_Hard_Spheres_get_diameter

    subroutine Hard_Spheres_set_data(this, snap_unit)

        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit

        character(len=4096) :: name
        integer :: num_particles
        integer :: snap_factor

        read(snap_unit, *) name, num_particles, snap_factor

        if (name /= this%name) write(error_unit, *) "Warning: names are different."
        if (num_particles /= this%num) error stop "Numbers of particles are different."
        if (snap_factor /= this%snap_factor) error stop "Snap factors are different."

    end subroutine Hard_Spheres_set_data

    subroutine Hard_Spheres_set_positions(this, Box_size, i_step, snap_unit, positions_set)

        class(Hard_Spheres), intent(inout) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit        
        logical, intent(out) :: positions_set

        integer :: i_particle
        real(DP), dimension(num_dimensions) :: position_i
        
        positions_set = .false.
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num
                read(snap_unit, *) position_i
                this%positions(:, i_particle) = modulo(position_i, Box_size)
            end do      
            positions_set = .true.      
        end if

    end subroutine Hard_Spheres_set_positions

    subroutine Dipolar_Hard_Spheres_set_moments(this, i_step, snap_unit, orientations_set)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit
        logical, intent(out) :: orientations_set

        integer :: i_particle
        
        orientations_set = .false.
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num
                read(snap_unit, *) this%moments(:, i_particle)
            end do
            orientations_set = .true.
        end if

    end subroutine Dipolar_Hard_Spheres_set_moments

    subroutine Hard_Spheres_set_test_num_particles(this, data_post_json, object_field)

        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: data_post_json
        character(len=*), intent(in) :: object_field

        call this%set_widom_num_particles(data_post_json, object_field)

        select type (this)
            type is (Dipolar_Hard_Spheres)
                call this%set_field_num_particles(data_post_json, object_field)
        end select

    end subroutine Hard_Spheres_set_test_num_particles

    subroutine Hard_Spheres_set_widom_num_particles(this, data_post_json, object_field)

        class(Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: data_post_json
        character(len=*), intent(in) :: object_field

        character(len=4096) :: data_field
        logical :: found

        data_field = "Chemical Potential."//object_field//".number of Widom particles"
        call data_post_json%get(data_field, this%widom_num_particles, found)
        call test_data_found(data_field, found)

    end subroutine Hard_Spheres_set_widom_num_particles
    
    subroutine Dipolar_Hard_Spheres_set_field_num_particles(this, data_post_json, object_field)

        class(Dipolar_Hard_Spheres), intent(inout) :: this
        type(json_file), intent(inout) :: data_post_json
        character(len=*), intent(in) :: object_field

        character(len=4096) :: data_field
        logical :: found

        data_field = "Local Field."//object_field//".External.number of test particles"
        call data_post_json%get(data_field, this%field_num_particles, found)
        call test_data_found(data_field, found)

    end subroutine Dipolar_Hard_Spheres_set_field_num_particles
    
    !> Write a report of the component in a file
    
    subroutine Hard_Spheres_write_report(this, num_particles, report_json)
    
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: num_particles
        type(json_value), pointer, intent(in) :: report_json

        type(json_value), pointer :: properties_json, snap_json
        real(DP) :: density, compacity, concentration

        call json_value_create(properties_json)
        call to_object(properties_json, "Properties")
        call json_value_add(report_json, properties_json)

        if (this%num > 0) then
            density = real(this%num + 1, DP) / this%volume ! cheating ? cf. Widom
            call json_value_add(properties_json, "density", density)
            compacity = 4._DP/3._DP*PI*(this%diameter/2._DP)**3 * density
            call json_value_add(properties_json, "compacity", compacity)
            concentration = real(this%num, DP) / real(num_particles, DP)
            call json_value_add(properties_json, "concentration", concentration)
        end if

        nullify(properties_json)

        call json_value_create(snap_json)
        call to_object(snap_json, "Snap")
        call json_value_add(report_json, snap_json)

        call json_value_add(snap_json, "factor", this%snap_factor)

        nullify(snap_json)
        
    end subroutine Hard_Spheres_write_report
    
    !> Take a snap shot of the configuration
    
    subroutine Hard_Spheres_write_data(this, snap_unit)
    
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: snap_unit
        
        write(snap_unit, *) this%name, this%num, this%snap_factor
        
    end subroutine Hard_Spheres_write_data
      
    subroutine Hard_Spheres_write_positions(this, i_step, snap_unit, double_precision)
        
        class(Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit
        logical, intent(in), optional :: double_precision
    
        integer :: i_particle
        character(len=len(double_precision_format)) :: output_format
        
        output_format = simple_precision_format
        
        if (present(double_precision)) then
            if (double_precision) then
                output_format = double_precision_format
            end if
        end if
        
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num
                write(snap_unit, output_format) this%positions(:, i_particle)
            end do
        end if

    end subroutine Hard_Spheres_write_positions
    
    subroutine Dipolar_Hard_Spheres_write_moments(this, i_step, snap_unit, double_precision)
        
        class(Dipolar_Hard_Spheres), intent(in) :: this
        integer, intent(in) :: i_step
        integer, intent(in) :: snap_unit
        logical, intent(in), optional :: double_precision
    
        integer :: i_particle
        character(len=len(double_precision_format)) :: output_format
        
        output_format = simple_precision_format
        
        if (present(double_precision)) then
            if (double_precision) then
                output_format = double_precision_format
            end if
        end if
        
        if (modulo(i_step, this%snap_factor) == 0) then
            do i_particle = 1, this%num
                write(snap_unit, output_format) this%moments(:, i_particle)
            end do
        end if

    end subroutine Dipolar_Hard_Spheres_write_moments
    
    !> Do an overlap test
    
    subroutine Hard_Spheres_test_overlap(this, Box_size)
    
        class(Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
    
        integer :: i_particle, j_particle
        real(DP), dimension(num_dimensions) :: position_i, position_j
        real(DP) :: distance_ij
    
        do j_particle = 1, this%num
            do i_particle = j_particle+1, this%num
            
                position_i(:) = this%positions(:, i_particle)
                position_j(:) = this%positions(:, j_particle)
                distance_ij = PBC_distance(Box_size, position_i, position_j)
                
                if (distance_ij < this%diameter) then
                    write(error_unit, *) this%name, " Overlap !", i_particle, j_particle
                    write(error_unit, *) " distance_ij = ", distance_ij
                    error stop
                end if
                    
            end do
        end do
        
        write(output_unit, *) this%name, ": Overlap test: OK !"
    
    end subroutine Hard_Spheres_test_overlap
    
    subroutine Between_Hard_Spheres_test_overlap(this, Box_size, type1, type2)
    
        class(Between_Hard_Spheres), intent(in) :: this
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: type1, type2
        
        integer :: type1_i_particle, type2_i_particle
        real(DP) :: distance_between
        real(DP), dimension(num_dimensions) :: type1_position, type2_position
        
        do type1_i_particle = 1, type1%get_num_particles()
            do type2_i_particle = 1, type2%get_num_particles()
                    
                type1_position(:) = type1%get_position(type1_i_particle)
                type2_position(:) = type2%get_position(type2_i_particle)
                distance_between = PBC_distance(Box_size, type1_position, type2_position)
                if (distance_between < this%diameter) then
                    write(error_unit, *) this%name, ": Overlap !", type1_i_particle, type2_i_particle
                    write(error_unit, *) " distance_between = ", distance_between
                    error stop
                end if

            end do
        end do

        write(output_unit, *) this%name, ": Overlap test: OK !"
    
    end subroutine Between_Hard_Spheres_test_overlap

end module class_hard_spheres
