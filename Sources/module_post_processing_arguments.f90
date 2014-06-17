module module_post_processing_arguments

use, intrinsic :: iso_fortran_env, only: error_unit

private
public argument_to_file

contains

    subroutine argument_to_file(arg_num, file_name, length)
        
        integer, intent(in) :: arg_num
        character(len=4096), intent(inout) :: file_name
        integer, intent(inout) :: length
        
        integer :: stat
        logical :: exist

        call get_command_argument(arg_num, file_name, length, stat)
        if (stat /= 0) error stop "error get_command_argument"
        inquire(file=file_name(1:length), exist=exist)
        if (.not. exist) then
            write(error_unit, *) "missing file: ", file_name(1:length)
            error stop
        end if
        
    end subroutine argument_to_file

end module module_post_processing_arguments
