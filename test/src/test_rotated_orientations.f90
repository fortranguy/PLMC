program test_rotated_orientations

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_orientations, only: Abstract_Orientations, Concrete_Orientations
use class_rotated_orientations, only: Abstract_Rotated_Orientations, Concrete_Rotated_Orientations

implicit none

    class(Abstract_Rotated_Orientations), allocatable :: rotated_orientations
    class(Abstract_Orientations), allocatable :: orientations
    class(Abstract_Particles_Number), allocatable :: particles_number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    integer :: num_steps, i_step, num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: orientation(:)
    real(DP) :: rotated_orientations_delta
    real(DP), dimension(num_dimensions) :: old_orientation, new_orientation
    integer, allocatable :: orientations_units(:)

    call json_initialize()
    data_filename = "rotated_orientations.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: particles_number)
    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call particles_number%set(num_particles)

    allocate(Concrete_Orientations :: orientations)
    call orientations%construct(particles_number)

    allocate(orientations_units(orientations%get_num()))
    do i_particle = 1, orientations%get_num()
        write(string_i, *) i_particle
        data_field = "Particles."//trim(adjustl(string_i))//".initial orientation"
        call input_data%get(data_field, orientation, data_found)
        call test_data_found(data_field, data_found)
        call orientations%set(i_particle, orientation)
        open(newunit=orientations_units(i_particle), recl=4096, &
             file="orientations_"//trim(adjustl(string_i))//".out", action="write")
    end do

    allocate(Concrete_Rotated_Orientations :: rotated_orientations)
    call rotated_orientations%construct(orientations)
    data_field = "Small Rotation.delta"
    call input_data%get(data_field, rotated_orientations_delta, data_found)
    call test_data_found(data_field, data_found)
    call rotated_orientations%set(rotated_orientations_delta)

    data_field = "Number of Steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)

    do i_step = 1, num_steps
        do i_particle = 1, orientations%get_num()
            old_orientation = orientations%get(i_particle)
            call orientations%set(i_particle, rotated_orientations%get(i_particle))
            new_orientation = orientations%get(i_particle)
            write(orientations_units(i_particle), *) i_step, old_orientation, &
                new_orientation - old_orientation
        end do
    end do

    do i_particle = orientations%get_num(), 1, -1
        close(orientations_units(i_particle))
    end do
    deallocate(orientations_units)

    call rotated_orientations%destroy()
    deallocate(rotated_orientations)
    call orientations%destroy()
    deallocate(orientations)
    deallocate(particles_number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_rotated_orientations
