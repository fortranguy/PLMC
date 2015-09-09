program test_field_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use types_field_parameters, only: Abstract_Field_Parameters, Constant_Field_Parameters
use class_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression
use procedures_errors, only: error_exit

implicit none

    class(Abstract_Field_Expression), allocatable :: field_expression
    class(Abstract_Field_Parameters), allocatable :: field_parameters
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field, field_name
    logical :: found
    real(DP), allocatable :: field(:), position(:)
    
    call json_initialize()
    data_filename = "field_expression.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    data_field = "Field.name"
    call input_data%get(data_field, field_name, found)
    call test_data_found(data_field, found)
    
    select case (field_name)
        case ("constant")
            allocate(Constant_Field_Parameters :: field_parameters)
            data_field = "Field.field"
            call input_data%get(data_field, field, found)
            call test_data_found(data_field, found)
            select type (field_parameters)
                type is (Constant_Field_Parameters)
                    field_parameters%field = field
            end select
            allocate(Constant_Field_Expression :: field_expression)
        case default
            call error_exit("Field not yet implemented?")
    end select
    
    call field_expression%set(field_parameters)
    data_field = "Field.position"
    call input_data%get(data_field, position, found)
    call test_data_found(data_field, found)
    write(output_unit, *) field_expression%get(position)
    
    deallocate(field_expression)
    deallocate(field_parameters)
    
    deallocate(field_name)
    deallocate(data_field)
    deallocate(data_filename)

end program test_field_expression
