module class_particles_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_coordinates, only: read_coordinates
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_number, only: Abstract_Particles_Number, &
    Concrete_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter, &
    Concrete_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, &
    Concrete_Particles_Moment_Norm, Null_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions, &
    Concrete_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations, &
    Concrete_Particles_Orientations, Null_Particles_Orientations
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade
use class_particles_total_moment, only: Particles_Total_Moment_Facade
use types_particles, only: Particles_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Particles_Factory
    private
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
    contains
        procedure :: allocate => Abstract_Particles_Factory_allocate
        procedure, private :: allocate_number =>  Abstract_Particles_Factory_allocate_number
        procedure, private :: allocate_diameter =>  Abstract_Particles_Factory_allocate_diameter
        procedure, private :: allocate_moment_norm => &
            Abstract_Particles_Factory_allocate_moment_norm
        procedure, private, nopass :: allocate_positions => &
            Abstract_Particles_Factory_allocate_positions
        procedure, private, nopass :: allocate_orientations => &
            Abstract_Particles_Factory_allocate_orientations
        procedure :: construct => Abstract_Particles_Factory_construct
        procedure, private :: construct_positions => &
            Abstract_Particles_Factory_construct_positions
        procedure, private :: construct_orientations => &
            Abstract_Particles_Factory_construct_orientations
        procedure :: destroy => Abstract_Particles_Factory_destroy
    end type Abstract_Particles_Factory

    type, extends(Abstract_Particles_Factory), public :: &
        Dipolar_Particles_Factory

    end type Dipolar_Particles_Factory

    type, extends(Abstract_Particles_Factory), public :: Apolar_Particles_Factory
    contains
        procedure, private :: allocate_moment_norm => &
            Apolar_Particles_Factory_allocate_moment_norm
        procedure, private, nopass :: allocate_orientations => &
            Apolar_Particles_Factory_allocate_orientations
    end type Apolar_Particles_Factory

contains

!implementation Abstract_Particles_Factory

    subroutine Abstract_Particles_Factory_allocate(this, particles, input_data, prefix)
        class(Abstract_Particles_Factory), intent(out) :: this
        type(Particles_Wrapper), intent(out) :: particles
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        this%input_data => input_data
        this%prefix = prefix
        call this%allocate_number(particles)
        call this%allocate_diameter(particles)
        call this%allocate_moment_norm(particles)
        call this%allocate_positions(particles)
        call this%allocate_orientations(particles)
    end subroutine Abstract_Particles_Factory_allocate

    subroutine Abstract_Particles_Factory_allocate_number(this, particles)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_particles

        data_field = this%prefix//".number"
        call this%input_data%get(data_field, num_particles, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Number :: particles%number)
        call particles%number%set(num_particles)
        deallocate(data_field)
    end subroutine Abstract_Particles_Factory_allocate_number

    subroutine Abstract_Particles_Factory_allocate_diameter(this, particles)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: diameter, diameter_min_factor

        data_field = this%prefix//".diameter"
        call this%input_data%get(data_field, diameter, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".minimum diameter factor"
        call this%input_data%get(data_field, diameter_min_factor, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Diameter :: particles%diameter)
        call particles%diameter%set(diameter, diameter_min_factor)
        deallocate(data_field)
    end subroutine Abstract_Particles_Factory_allocate_diameter

    subroutine Abstract_Particles_Factory_allocate_moment_norm(this, particles)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        data_field = this%prefix//".moment norm"
        call this%input_data%get(data_field, moment_norm, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Moment_Norm :: particles%moment_norm)
        call particles%moment_norm%set(moment_norm)
        deallocate(data_field)
    end subroutine Abstract_Particles_Factory_allocate_moment_norm

    subroutine Abstract_Particles_Factory_allocate_positions(particles)
        type(Particles_Wrapper), intent(inout) :: particles

        allocate(Concrete_Particles_Positions :: particles%positions)
    end subroutine Abstract_Particles_Factory_allocate_positions

    subroutine Abstract_Particles_Factory_allocate_orientations(particles)
        type(Particles_Wrapper), intent(inout) :: particles

        allocate(Concrete_Particles_Orientations :: particles%orientations)
    end subroutine Abstract_Particles_Factory_allocate_orientations

    subroutine Abstract_Particles_Factory_construct(this, particles, periodic_box)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call this%construct_positions(particles, periodic_box)
        call this%construct_orientations(particles)
        call particles%dipolar_moments%construct(particles%moment_norm, particles%orientations)
        call particles%total_moment%construct(particles%dipolar_moments)
    end subroutine Abstract_Particles_Factory_construct

    subroutine Abstract_Particles_Factory_construct_positions(this, particles, periodic_box)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field, filename
        logical :: data_found
        real(DP), allocatable :: positions(:, :)
        integer :: i_particle

        call particles%positions%construct(periodic_box, particles%number)
        data_field = this%prefix//".initial positions"
        call this%input_data%get(data_field, filename, data_found)
        call test_data_found(data_field, data_found)
        call read_coordinates(positions, particles%positions%get_num(), filename)
        do i_particle = 1, particles%positions%get_num()
            call particles%positions%set(i_particle, positions(:, i_particle))
        end do
        if (allocated(positions)) deallocate(positions)
        deallocate(filename)
        deallocate(data_field)
    end subroutine Abstract_Particles_Factory_construct_positions

    subroutine Abstract_Particles_Factory_construct_orientations(this, particles)
        class(Abstract_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        character(len=:), allocatable :: data_field, filename
        logical :: data_found
        real(DP), allocatable :: orientations(:, :)
        integer :: i_particle

        call particles%orientations%construct(particles%number)
        data_field = this%prefix//".initial orientations"
        call this%input_data%get(data_field, filename, data_found)
        call test_data_found(data_field, data_found)
        call read_coordinates(orientations, particles%orientations%get_num(), filename)
        do i_particle = 1, particles%orientations%get_num()
            call particles%orientations%set(i_particle, orientations(:, i_particle))
        end do
        if (allocated(orientations)) deallocate(orientations)
        deallocate(filename)
        deallocate(data_field)
    end subroutine Abstract_Particles_Factory_construct_orientations

    subroutine Abstract_Particles_Factory_destroy(this, particles)
        class(Abstract_Particles_Factory), intent(inout) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        call particles%orientations%destroy()
        if (allocated(particles%orientations)) deallocate(particles%orientations)
        call particles%positions%destroy()
        if (allocated(particles%positions)) deallocate(particles%positions)
        if (allocated(particles%moment_norm)) deallocate(particles%moment_norm)
        if (allocated(particles%diameter)) deallocate(particles%diameter)
        if (allocated(particles%number)) deallocate(particles%number)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Abstract_Particles_Factory_destroy

!end implementation Abstract_Particles_Factory

!implementation Apolar_Particles_Factory

    subroutine Apolar_Particles_Factory_allocate_moment_norm(this, particles)
        class(Apolar_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        allocate(Null_Particles_Moment_Norm :: particles%moment_norm)
    end subroutine Apolar_Particles_Factory_allocate_moment_norm

    subroutine Apolar_Particles_Factory_allocate_orientations(particles)
        type(Particles_Wrapper), intent(inout) :: particles

        allocate(Null_Particles_Orientations :: particles%orientations)
    end subroutine Apolar_Particles_Factory_allocate_orientations

!end implementation Apolar_Particles_Factory

end module class_particles_factory
