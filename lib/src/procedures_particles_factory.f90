module procedures_particles_factory

use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_dipolar_moments, only: Particles_Dipolar_Moments_Facade
use class_particles_total_moment, only: Particles_Total_Moment_Facade
use types_particles, only: Concrete_Particles

implicit none

private
public :: particles_construct, particles_destroy

contains

    subroutine particles_construct(particles, input_data, prefix)
        type(Concrete_Particles), intent(out) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call set_number(particles, input_data, prefix)
        !particles%diameter => diameter
        !particles%moment_norm => moment_norm
        !particles%positions => positions
        !particles%orientations => orientations
    end subroutine particles_construct

    subroutine particles_destroy(particles)
        type(Concrete_Particles), intent(inout) :: particles

        !particles%orientations
        !particles%positions
        !particles%moment_norm
        !particles%diameter
        if (allocated(particles%number)) deallocate(particles%number)
    end subroutine particles_destroy

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

end module procedures_particles_factory
