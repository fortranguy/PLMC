program test_sphere

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit

use module_data, only: test_file_exists, test_data_found
use class_sphere, only : Abstract_Sphere, Concrete_Sphere
use json_module, only: json_file, json_initialize

implicit none

    class(Abstract_Sphere), allocatable :: sphere
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP) :: diameter

    call json_initialize()

    data_filename = "sphere.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Sphere :: sphere)

    data_field = "Sphere.diameter"
    call input_data%get(data_field, diameter, found)
    call test_data_found(data_field, found)

    call sphere%set_diameter(diameter)
    write(output_unit, *) "diameter =", sphere%get_diameter()

    deallocate(sphere)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_sphere