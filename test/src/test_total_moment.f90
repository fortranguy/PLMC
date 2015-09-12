program test_total_moment

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_orientation, only: random_orientation
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
    logical :: data_found, write_orientation
    integer ::num_particles, i_particle
    real(DP) :: moment_norm
    integer ::dipolar_moments_unit

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
    data_field = "Particles.moment norm"
    call input_data%get(data_field, moment_norm, data_found)
    call test_data_found(data_field, data_found)
    do i_particle = 1, moments_norm%get_num()
        call moments_norm%set(i_particle, moment_norm)
    end do

    allocate(Concrete_Orientations :: orientations)
    call orientations%construct(particles_number)
    do i_particle = 1, orientations%get_num()
        call orientations%set(i_particle, random_orientation())
    end do

    call dipolar_moments%construct(moments_norm, orientations)
    data_field = "Write moments"
    call input_data%get(data_field, write_orientation, data_found)
    call test_data_found(data_field, data_found)
    if (write_orientation) then
        open(newunit=dipolar_moments_unit, recl=4096, file="dipolar_moments.out", action="write")
        do i_particle = 1, dipolar_moments%get_num()
            write(dipolar_moments_unit, *) dipolar_moments%get(i_particle)
        end do
        close(dipolar_moments_unit)
    end if

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
