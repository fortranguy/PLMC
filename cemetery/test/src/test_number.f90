program test_number

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists, check_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number

implicit none

    class(Abstract_Particles_Number), allocatable :: number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: num_particles

    call json_initialize()

    data_filename = "number.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: number)

    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, found)
    call check_data_found(data_field, found)

    call number%set(num_particles)
    write(output_unit, *) "number =", number%get()

    deallocate(number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_number
