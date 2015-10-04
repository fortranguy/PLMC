module procedures_command_arguments

use procedures_errors, only: error_exit
use module_data, only: test_file_exists

implicit none
private
public :: set_filename_from_argument

contains

    subroutine set_filename_from_argument(filename)
        character(len=:), allocatable, intent(out) :: filename

        if (command_argument_count() /= 1) then
            call error_exit("Please provide a json input file as the only argument.")
        end if
        call argument_to_file(1, filename)
    end subroutine set_filename_from_argument

    subroutine argument_to_file(i_argument, filename_i)
        integer, intent(in) :: i_argument
        character(len=:), allocatable, intent(out) :: filename_i

        character(len=1024) :: filename
        integer :: length, argument_stat

        call get_command_argument(i_argument, filename, length, argument_stat)
        if (argument_stat /= 0) then
            call error_exit("argument_to_file: error")
        end if
        call test_file_exists(filename(1:length))
        filename_i = filename(1:length)
    end subroutine argument_to_file

end module procedures_command_arguments
