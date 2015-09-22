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
use class_particles_chemical_potential, only : Abstract_Particles_Chemical_Potential, &
    Concrete_Particles_Chemical_Potential
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade
use class_particles_total_moment, only: Particles_Total_Moment_Facade
use types_particles, only: Concrete_Particles_Parameters, Particles_Wrapper

implicit none

private

    type, abstract, public :: Concrete_Particles_Factory
    private
        type(Concrete_Particles_Parameters) :: parameters
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
    contains
        procedure :: allocate => Concrete_Particles_Factory_allocate
        procedure, private :: allocate_number =>  Concrete_Particles_Factory_allocate_number
        procedure, private :: allocate_diameter =>  Concrete_Particles_Factory_allocate_diameter
        procedure, private :: allocate_moment_norm => &
            Concrete_Particles_Factory_allocate_moment_norm
        procedure, private :: allocate_positions => Concrete_Particles_Factory_allocate_positions
        procedure, private :: allocate_orientations => &
            Concrete_Particles_Factory_allocate_orientations
        procedure, private :: allocate_chemical_potential => &
            Concrete_Particles_Factory_allocate_chemical_potential
        procedure :: construct => Concrete_Particles_Factory_construct
        procedure, private :: construct_positions => &
            Concrete_Particles_Factory_construct_positions
        procedure, private :: construct_orientations => &
            Concrete_Particles_Factory_construct_orientations
        procedure :: destroy => Concrete_Particles_Factory_destroy
    end type Concrete_Particles_Factory

contains

!implementation Concrete_Particles_Factory

    subroutine Concrete_Particles_Factory_allocate(this, parameters, particles, input_data, prefix)
        class(Concrete_Particles_Factory), intent(out) :: this
        type(Concrete_Particles_Parameters), intent(in) :: parameters
        type(Particles_Wrapper), intent(out) :: particles
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        this%parameters%exist = parameters%exist
        this%parameters%are_dipolar = parameters%are_dipolar
        this%parameters%can_exchange = parameters%can_exchange
        this%input_data => input_data
        this%prefix = prefix
        call this%allocate_number(particles%number)
        call this%allocate_diameter(particles%diameter)
        call this%allocate_moment_norm(particles%moment_norm)
        call this%allocate_positions(particles%positions)
        call this%allocate_orientations(particles%orientations)
        call this%allocate_chemical_potential(particles%chemical_potential)
    end subroutine Concrete_Particles_Factory_allocate

    subroutine Concrete_Particles_Factory_allocate_number(this, particles_number)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Number), allocatable, intent(out) :: particles_number

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_particles

        data_field = this%prefix//".number"
        call this%input_data%get(data_field, num_particles, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Number :: particles_number)
        call particles_number%set(num_particles)
        deallocate(data_field)
    end subroutine Concrete_Particles_Factory_allocate_number

    subroutine Concrete_Particles_Factory_allocate_diameter(this, particles_diameter)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Diameter), allocatable, intent(out) :: particles_diameter

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: diameter, diameter_min_factor

        data_field = this%prefix//".diameter"
        call this%input_data%get(data_field, diameter, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".minimum diameter factor"
        call this%input_data%get(data_field, diameter_min_factor, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Diameter :: particles_diameter)
        call particles_diameter%set(diameter, diameter_min_factor)
        deallocate(data_field)
    end subroutine Concrete_Particles_Factory_allocate_diameter

    subroutine Concrete_Particles_Factory_allocate_moment_norm(this, particles_moment_norm)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Moment_Norm), allocatable, intent(out) :: particles_moment_norm

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        data_field = this%prefix//".moment norm"
        call this%input_data%get(data_field, moment_norm, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Moment_Norm :: particles_moment_norm)
        call particles_moment_norm%set(moment_norm)
        deallocate(data_field)
    end subroutine Concrete_Particles_Factory_allocate_moment_norm

    subroutine Concrete_Particles_Factory_allocate_positions(this, particles_positions)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Positions), allocatable, intent(out) :: particles_positions

        allocate(Concrete_Particles_Positions :: particles_positions)
    end subroutine Concrete_Particles_Factory_allocate_positions

    subroutine Concrete_Particles_Factory_allocate_orientations(this, particles_orientations)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Orientations), allocatable, intent(out) :: particles_orientations

        allocate(Concrete_Particles_Orientations :: particles_orientations)
    end subroutine Concrete_Particles_Factory_allocate_orientations

    subroutine Concrete_Particles_Factory_allocate_chemical_potential(this, &
        particles_chemical_potential)
        class(Concrete_Particles_Factory), intent(in) :: this
        class(Abstract_Particles_Chemical_Potential), allocatable, intent(out) :: &
            particles_chemical_potential

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, excess

        data_field = this%prefix//"Chemical Potential.density"
        call this%input_data%get(data_field, density, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//"Chemical Potential.excess"
        call this%input_data%get(data_field, excess, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Chemical_Potential :: particles_chemical_potential)
        call particles_chemical_potential%set(density, excess)
        deallocate(data_field)
    end subroutine Concrete_Particles_Factory_allocate_chemical_potential

    subroutine Concrete_Particles_Factory_construct(this, particles, periodic_box)
        class(Concrete_Particles_Factory), intent(in) :: this
        type(Particles_Wrapper), intent(inout) :: particles
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call this%construct_positions(particles, periodic_box)
        call this%construct_orientations(particles)
        call particles%dipolar_moments%construct(particles%moment_norm, particles%orientations)
        call particles%total_moment%construct(particles%dipolar_moments)
    end subroutine Concrete_Particles_Factory_construct

    subroutine Concrete_Particles_Factory_construct_positions(this, particles, periodic_box)
        class(Concrete_Particles_Factory), intent(in) :: this
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
    end subroutine Concrete_Particles_Factory_construct_positions

    subroutine Concrete_Particles_Factory_construct_orientations(this, particles)
        class(Concrete_Particles_Factory), intent(in) :: this
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
    end subroutine Concrete_Particles_Factory_construct_orientations

    subroutine Concrete_Particles_Factory_destroy(this, particles)
        class(Concrete_Particles_Factory), intent(inout) :: this
        type(Particles_Wrapper), intent(inout) :: particles

        if (allocated(particles%chemical_potential)) deallocate(particles%chemical_potential)
        call particles%orientations%destroy()
        if (allocated(particles%orientations)) deallocate(particles%orientations)
        call particles%positions%destroy()
        if (allocated(particles%positions)) deallocate(particles%positions)
        if (allocated(particles%moment_norm)) deallocate(particles%moment_norm)
        if (allocated(particles%diameter)) deallocate(particles%diameter)
        if (allocated(particles%number)) deallocate(particles%number)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Concrete_Particles_Factory_destroy

!end implementation Concrete_Particles_Factory

end module class_particles_factory
