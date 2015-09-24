program test_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_field_expression, only: Abstract_Field_Expression
use procedures_box_factory, only: allocate_and_set_field_expression

implicit none

    class(Abstract_Field_Expression), allocatable :: field_expression
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: position(:)

    call json_initialize()
    data_filename = "field_expression.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call allocate_and_set_field_expression(field_expression, input_data, "Test Field Expression")

    data_field = "Test Field Expression.Particle.position"
    call input_data%get(data_field, position, found)
    call test_data_found(data_field, found)
    deallocate(data_field)
    write(output_unit, *) field_expression%get(position)
    deallocate(position)

    deallocate(field_expression)

end program test_field_expression
