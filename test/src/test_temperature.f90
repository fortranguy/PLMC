program test_temperature

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists, check_data_found
use class_temperature, only: Abstract_Temperature, Concrete_Temperature

implicit none

    class(Abstract_Temperature), allocatable :: temperature
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP) :: temperature_value
    
    call json_initialize()
    data_filename = "temperature.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    data_field = "Temperature.value"
    allocate(Concrete_Temperature :: temperature)
    call input_data%get(data_field, temperature_value, found)
    call check_data_found(data_field, found)
    call temperature%set(temperature_value)
    write(output_unit, *) "temperature =", temperature%get()
    
    deallocate(temperature)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_temperature
