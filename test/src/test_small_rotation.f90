program test_small_rotation

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_orientations, only: Abstract_Orientations, Concrete_Orientations
use class_small_rotation, only: Abstract_Small_Rotation, Concrete_Small_Rotation

implicit none

    class(Abstract_Small_Rotation), allocatable :: small_rotation
    class(Abstract_Orientations), allocatable :: orientations
    class(Abstract_Particles_Number), allocatable :: particles_number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    integer :: num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: orientation(:)
    real(DP) :: small_rotation_delta

    integer :: orientations_unit

    call json_initialize()
    data_filename = "small_rotation.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: particles_number)
    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call particles_number%set(num_particles)

    allocate(Concrete_Orientations :: orientations)
    call orientations%construct(particles_number)

    do i_particle = 1, particles_number%get()
        write(string_i, *) i_particle
        data_field = "Particles."//trim(adjustl(string_i))//".orientation"
        call input_data%get(data_field, orientation, data_found)
        call test_data_found(data_field, data_found)
        call orientations%set(i_particle, orientation)
    end do

    allocate(Concrete_Small_Rotation :: small_rotation)
    call small_rotation%construct(orientations)
    data_field = "Small Rotation.delta"
    call input_data%get(data_field, small_rotation_delta, data_found)
    call test_data_found(data_field, data_found)
    call small_rotation%set(small_rotation_delta)

    open(newunit=orientations_unit, recl=4096, file="rotated_orientations.out", action="write")
    do i_particle = 1, particles_number%get()
        write(orientations_unit, *) orientations%get(i_particle), &
            small_rotation%get(i_particle) - orientations%get(i_particle)
        write(orientations_unit, *)
        write(orientations_unit, *)
    end do
    close(orientations_unit)

    call small_rotation%destroy()
    deallocate(small_rotation)
    call orientations%destroy()
    deallocate(orientations)
    deallocate(particles_number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_small_rotation
