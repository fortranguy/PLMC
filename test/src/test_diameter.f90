program test_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_diameter, only: Abstract_Diameter, Concrete_Diameter

implicit none

    class(Abstract_Diameter), allocatable :: diameter
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    real(DP) :: diameter_value, diameter_min_factor

    call json_initialize()

    data_filename = "diameter.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    allocate(Concrete_Diameter :: diameter)
    data_field = "Diameter.value"
    call input_data%get(data_field, diameter_value, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Diameter.minimum factor"
    call input_data%get(data_field, diameter_min_factor, data_found)
    call test_data_found(data_field, data_found)
    call diameter%set(diameter_value, diameter_min_factor)
    deallocate(data_field)

    write(output_unit, *) "diameter", diameter%get()
    write(output_unit, *) "minimum diameter", diameter%get_min()

    deallocate(diameter)
    call input_data%destroy()

end program test_diameter