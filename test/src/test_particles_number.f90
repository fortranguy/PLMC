program test_particles_number

use, intrinsic :: iso_fortran_env, only: output_unit

use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Particles_Number), allocatable :: particles_num
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: num

    call json_initialize()

    data_filename = "particles_number.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Particles_Number :: particles_num)

    data_field = "Particles Number.number"
    call input_data%get(data_field, num, found)
    call test_data_found(data_field, found)

    call particles_num%set_num(num)
    write(output_unit, *) "number =", particles_num%get_num()

    deallocate(particles_num)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_particles_number