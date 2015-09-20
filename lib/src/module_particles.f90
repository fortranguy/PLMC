module module_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade
use class_particles_total_moment, only: Particles_Total_Moment_Facade

implicit none

private
public :: Concrete_Particles_construct, Concrete_Particles_destroy

    type, public :: Concrete_Particle
        logical :: same_type = .false.
        integer :: i = 0
        real(DP) :: diameter = 0._DP
        real(DP) :: min_diameter = 0._DP
        real(DP) :: moment_norm = 0._DP
        real(DP) :: position(num_dimensions) = 0._DP
        real(DP) :: orientation(num_dimensions) = 0._DP
    end type Concrete_Particle

    type, public :: Concrete_Particles
        class(Abstract_Particles_Number), allocatable :: number
        class(Abstract_Particles_Diameter), allocatable :: diameter
        class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
        class(Abstract_Particles_Positions), allocatable :: positions
        class(Abstract_Particles_Orientations), allocatable :: orientations
        class(Particles_Dipolar_Moments_Facade), allocatable :: dipolar_moments
        class(Particles_Total_Moment_Facade), allocatable :: total_moment
    end type Concrete_Particles

contains

    subroutine Concrete_Particles_construct(particles, input_data, prefix)
        type(Concrete_Particles), intent(out) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        !particles%number => number
        call set_number(particles, input_data, prefix)
        !particles%diameter => diameter
        !particles%moment_norm => moment_norm
        !particles%positions => positions
        !particles%orientations => orientations
    end subroutine Concrete_Particles_construct

    subroutine Concrete_Particles_destroy(particles)
        type(Concrete_Particles), intent(inout) :: particles

        !particles%orientations
        !particles%positions
        !particles%moment_norm
        !particles%diameter
        !particles%number
    end subroutine Concrete_Particles_destroy

    subroutine set_number(particles, input_data, prefix)
        type(Concrete_Particles), intent(inout) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_particles

        data_field = prefix//".number"
        call input_data%get(data_field, num_particles, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Particles_Number :: particles%number)
        call particles%number%set(num_particles)
        deallocate(data_field)
    end subroutine set_number

end module module_particles
