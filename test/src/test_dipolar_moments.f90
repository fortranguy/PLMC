module procedures_dipolar_moments_print

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments

implicit none

private
public set_objects, print_dipolar_moments

contains

    subroutine set_objects(input_data, moment_norm_name, moment_norm, &
                           orientations_name, orientations)
        type(json_file), intent(inout) :: input_data
        class(Abstract_Particles_Moment_Norm), intent(inout) :: moment_norm
        class(Abstract_Particles_Orientations), intent(inout) :: orientations
        character(len=*), intent(in) :: moment_norm_name, orientations_name

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: i_particle
        character(len=1024) :: string_i
        real(DP) :: moment_norm_value
        real(DP), allocatable :: orientation(:)

        write(output_unit, *) moment_norm_name//" + "//orientations_name

        data_field = "Dipoles.moment norm"
        call input_data%get(data_field, moment_norm_value, data_found)
        call test_data_found(data_field, data_found)
        call moment_norm%set(moment_norm_value)

        do i_particle = 1, orientations%get_num()
            write(string_i, *) i_particle
            data_field = "Dipoles.orientation "//trim(adjustl(string_i))
            call input_data%get(data_field, orientation, data_found)
            call test_data_found(data_field, data_found)
            call orientations%set(i_particle, orientation)
            deallocate(orientation)
        end do
    end subroutine set_objects

    subroutine print_dipolar_moments(dipolar_moments)
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: dipolar_moments

        integer :: i_particle

        do i_particle = 1, dipolar_moments%get_num()
            write(output_unit, *) dipolar_moments%get(i_particle)
        end do
    end subroutine print_dipolar_moments

end module procedures_dipolar_moments_print

program test_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, &
    Null_Particles_Moment_Norm, Concrete_Particles_Moment_Norm
use class_particles_orientations, only: Abstract_Particles_Orientations, &
    Null_Particles_Orientations, Concrete_Particles_Orientations
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments, &
    Concrete_Particles_Dipolar_Moments
use procedures_dipolar_moments_print, only: set_objects, print_dipolar_moments

implicit none

    class(Abstract_Particles_Dipolar_Moments), allocatable :: dipolar_moments
    class(Abstract_Particles_Orientations), allocatable :: orientations
    class(Abstract_Particles_Moment_Norm), allocatable :: moment_norm
    class(Abstract_Particles_Number), allocatable :: number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    integer :: num_particles

    call json_initialize()
    data_filename = "dipolar_moments.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: number)
    data_field = "Dipoles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call number%set(num_particles)

    allocate(Null_Particles_Moment_Norm :: moment_norm)
    allocate(Null_Particles_Orientations :: orientations)
    call orientations%construct(number)
    call set_objects(input_data, "Null_Particles_Moment_Norm", moment_norm, &
                     "Null_Particles_Orientations", orientations)
    allocate(Concrete_Particles_Dipolar_Moments :: dipolar_moments)
    call dipolar_moments%construct(moment_norm, orientations)
    call print_dipolar_moments(dipolar_moments)
    call dipolar_moments%destroy()
    deallocate(dipolar_moments)
    call orientations%destroy()
    deallocate(orientations)
    deallocate(moment_norm)

    allocate(Concrete_Particles_Moment_Norm :: moment_norm)
    allocate(Concrete_Particles_Orientations :: orientations)
    call orientations%construct(number)
    call set_objects(input_data, "Uniform_Moments", moment_norm, &
                     "Concrete_Particles_Orientations", orientations)
    call dipolar_moments%construct(moment_norm, orientations)
    call print_dipolar_moments(dipolar_moments)
    call dipolar_moments%destroy()
    call orientations%destroy()
    deallocate(orientations)
    deallocate(moment_norm)

    deallocate(number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_dipolar_moments
