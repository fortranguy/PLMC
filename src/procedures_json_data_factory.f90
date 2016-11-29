module procedures_json_data_factory

use, intrinsic :: iso_fortran_env, only: error_unit
use json_module, only: json_file, json_value
use procedures_command_arguments, only: create_filename_from_argument

implicit none

private
public :: create, destroy

contains

    subroutine create(input_data, i_argument)
        type(json_file), intent(out) :: input_data
        integer, intent(in) :: i_argument

        character(len=:), allocatable :: data_filename

        call input_data%initialize()
        if (input_data%failed()) call input_data%print_error_message(error_unit)
        call create_filename_from_argument(data_filename, i_argument)
        call input_data%load_file(data_filename)
        if (input_data%failed()) call input_data%print_error_message(error_unit)
    end subroutine create

     subroutine destroy(input_data)
        type(json_file), intent(inout) :: input_data

        call input_data%destroy()
        if (input_data%failed()) call input_data%print_error_message(error_unit)
    end subroutine destroy

end module procedures_json_data_factory
