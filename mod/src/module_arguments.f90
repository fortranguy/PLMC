module module_arguments

use, intrinsic :: iso_fortran_env, only: error_unit

implicit none
private
public arg_to_file

contains

    subroutine arg_to_file(arg_num, filename, length)
        
        integer, intent(in) :: arg_num
        character(len=4096), intent(inout) :: filename
        integer, intent(inout) :: length
        
        integer :: stat
        logical :: exist

        call get_command_argument(arg_num, filename, length, stat)
        if (stat /= 0) error stop "error get_command_argument"
        inquire(file=filename(1:length), exist=exist)
        if (.not. exist) then
            write(error_unit, *) "missing file: ", filename(1:length)
            error stop
        end if
        
    end subroutine arg_to_file

end module module_arguments
