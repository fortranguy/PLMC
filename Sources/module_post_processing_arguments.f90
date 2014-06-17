module module_post_processing_arguments

use, intrinsic :: iso_fortran_env, only: error_unit

private
public argument_to_file

contains

    subroutine argument_to_file(number, file, length)
        
        integer, intent(in) :: number
        character(len=4096), intent(inout) :: file
        integer, intent(inout) :: length
        
        integer :: stat
        logical :: exist

        call get_command_argument(number, file, length, stat)
        if (stat /= 0) error stop "error get_command_argument"
        inquire(file=file(1:length), exist=exist)
        if (.not. exist) then
            write(error_unit, *) "missing file: ", file(1:length)
            error stop
        end if
        
    end subroutine argument_to_file

end module module_post_processing_arguments
