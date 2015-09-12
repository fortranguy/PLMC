program test_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_moments_norm, only: Abstract_Moments_Norm, Uniform_Moments_Norm
use class_orientations, only: Abstract_Orientations, Concrete_Orientations
use class_dipolar_moments, only: Dipolar_Moments_Facade
use class_total_moment, only: Abstract_Total_Moment, Concrete_Total_Moment

implicit none

    class(Abstract_Total_Moment), allocatable :: total_moment
    type(Dipolar_Moments_Facade) :: dipolar_moments
    class(Abstract_Orientations), allocatable :: orientations
    class(Abstract_Moments_Norm), allocatable :: moments_norm
    class(Abstract_Particles_Number), allocatable :: particles_number

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    integer ::num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: orientation(:)
    real(DP) :: moment_norm

    call json_initialize()
    data_filename = "total_moment.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: particles_number)
    data_field = "Number of Particles"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call particles_number%set(num_particles)

    allocate(Uniform_Moments_Norm :: moments_norm)
    call moments_norm%construct(particles_number)
    do i_particle = 1, moments_norm%get_num()
        write(string_i, *) i_particle
        data_field = "Particle "//trim(adjustl(string_i))//".moment norm"
        call input_data%get(data_field, moment_norm, data_found)
        call test_data_found(data_field, data_found)
        call moments_norm%set(i_particle, moment_norm)
    end do

    allocate(Concrete_Orientations :: orientations)
    call orientations%construct(particles_number)
    do i_particle = 1, orientations%get_num()
        write(string_i, *) i_particle
        data_field = "Particle "//trim(adjustl(string_i))//".orientation"
        call input_data%get(data_field, orientation, data_found)
        call test_data_found(data_field, data_found)
        call orientations%set(i_particle, orientation)
        deallocate(orientation)
    end do

    call dipolar_moments%construct(moments_norm, orientations)
    allocate(Concrete_Total_Moment :: total_moment)
    call total_moment%construct(dipolar_moments)

    write(output_unit, *) "Total Moment =", total_moment%get()

    call total_moment%destroy()
    deallocate(total_moment)
    call dipolar_moments%destroy()
    call orientations%destroy()
    deallocate(orientations)
    call moments_norm%destroy()
    deallocate(moments_norm)
    deallocate(particles_number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_total_moment
