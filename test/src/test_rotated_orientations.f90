program test_rotated_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_particles_orientations, only: Abstract_Particles_Orientations, &
    Concrete_Particles_Orientations
use class_rotated_particles_orientations, only: Abstract_Rotated_Particles_Orientations, &
    Concrete_Rotated_Particles_Orientations

implicit none

    class(Abstract_Rotated_Particles_Orientations), allocatable :: rotated_orientations
    class(Abstract_Particles_Orientations), allocatable :: orientations
    class(Abstract_Particles_Number), allocatable :: number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    integer :: num_steps, i_step, num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: orientation(:)
    real(DP) :: rotated_orientations_delta, adaptation_factor
    real(DP), dimension(num_dimensions) :: old_orientation, new_orientation
    integer, allocatable :: orientations_big_units(:), orientations_small_units(:)

    call json_initialize()
    data_filename = "rotated_orientations.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: number)
    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call number%set(num_particles)

    allocate(Concrete_Particles_Orientations :: orientations)
    call orientations%construct(number)

    allocate(orientations_big_units(orientations%get_num()))
    allocate(orientations_small_units(orientations%get_num()))
    do i_particle = 1, orientations%get_num()
        write(string_i, *) i_particle
        data_field = "Particles."//trim(adjustl(string_i))//".initial orientation"
        call input_data%get(data_field, orientation, data_found)
        call test_data_found(data_field, data_found)
        call orientations%set(i_particle, orientation)
        open(newunit=orientations_big_units(i_particle), recl=4096, &
             file="orientations_big_"//trim(adjustl(string_i))//".out", action="write")
        open(newunit=orientations_small_units(i_particle), recl=4096, &
             file="orientations_small_"//trim(adjustl(string_i))//".out", action="write")
    end do

    allocate(Concrete_Rotated_Particles_Orientations :: rotated_orientations)
    data_field = "Small Rotation.delta"
    call input_data%get(data_field, rotated_orientations_delta, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Small Rotation.adaptation factor"
    call input_data%get(data_field, adaptation_factor, data_found)
    call test_data_found(data_field, data_found)
    call rotated_orientations%construct(orientations, rotated_orientations_delta, adaptation_factor)

    data_field = "Number of Steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)

    do i_step = 1, num_steps
        call rotated_orientations%increase()
        do i_particle = 1, orientations%get_num()
            old_orientation = orientations%get(i_particle)
            call orientations%set(i_particle, rotated_orientations%get(i_particle))
            new_orientation = orientations%get(i_particle)
            write(orientations_big_units(i_particle), *) i_step, old_orientation, &
                new_orientation - old_orientation
        end do
    end do

    do i_particle = orientations%get_num(), 1, -1
        close(orientations_big_units(i_particle))
    end do
    deallocate(orientations_big_units)

    do i_step = 1, num_steps
        call rotated_orientations%decrease()
        do i_particle = 1, orientations%get_num()
            old_orientation = orientations%get(i_particle)
            call orientations%set(i_particle, rotated_orientations%get(i_particle))
            new_orientation = orientations%get(i_particle)
            write(orientations_small_units(i_particle), *) i_step, old_orientation, &
                new_orientation - old_orientation
        end do
    end do

    do i_particle = orientations%get_num(), 1, -1
        close(orientations_small_units(i_particle))
    end do
    deallocate(orientations_small_units)

    call rotated_orientations%destroy()
    deallocate(rotated_orientations)
    call orientations%destroy()
    deallocate(orientations)
    deallocate(number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_rotated_orientations
