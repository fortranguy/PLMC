program test_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists, check_data_found
use procedures_errors, only: error_exit
use procedures_property_inquirers, only: apply_external_field
use class_field_expression, only: Abstract_Field_Expression
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy

implicit none

    class(Abstract_Field_Expression), allocatable :: field_expression
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    real(DP), allocatable :: position(:)
    logical :: field_applied

    call json_initialize()
    data_filename = "field_expression.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    field_applied = apply_external_field(input_data, "Environment.")
    call environment_factory_create(field_expression, field_applied, input_data, &
        "Test Field Expression")

    data_field = "Test Field Expression.Particle.position"
    call input_data%get(data_field, position, found)
    call check_data_found(data_field, found)
    deallocate(data_field)
    write(output_unit, *) field_expression%get(position)
    deallocate(position)

    call environment_factory_destroy(field_expression)

end program test_field_expression
