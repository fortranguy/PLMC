module procedures_json_data_factory

use, intrinsic :: iso_fortran_env, only: error_unit
use json_module, only: json_core, json_file, json_value
use procedures_command_arguments, only: create_filename_from_argument

implicit none

private
public :: create_input, create_output, destroy_input, destroy_output

contains

    subroutine create_input(input_data, i_argument)
        type(json_file), intent(out) :: input_data
        integer, intent(in) :: i_argument

        character(len=:), allocatable :: data_filename

        call input_data%initialize()
        if (input_data%failed()) call input_data%print_error_message(error_unit)
        call create_filename_from_argument(data_filename, i_argument)
        call input_data%load_file(data_filename)
        if (input_data%failed()) call input_data%print_error_message(error_unit)
    end subroutine create_input

     subroutine destroy_input(input_data)
        type(json_file), intent(inout) :: input_data

        call input_data%destroy()
    end subroutine destroy_input

    subroutine create_output(json, output_data)
        type(json_core), intent(out) :: json
        type(json_value), pointer, intent(out) :: output_data

        call json%initialize()
        call json%create_object(output_data, "")
    end subroutine create_output

    subroutine destroy_output(json, output_data)
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(inout) :: output_data

        call json%destroy(output_data)
        if (json%failed()) call json%print_error_message(error_unit)
    end subroutine destroy_output

end module procedures_json_data_factory
