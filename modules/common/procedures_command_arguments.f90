module procedures_command_arguments

use data_constants, only: max_word_length
use procedures_errors, only: error_exit
use procedures_checks, only: check_file_exists

implicit none

private
public :: create_filename_from_argument

contains

    subroutine create_filename_from_argument(filename_i, i_argument)
        character(len=:), allocatable, intent(out) :: filename_i
        integer, intent(in) :: i_argument

        character(len=max_word_length) :: filename
        integer :: length, argument_stat

        call get_command_argument(i_argument, filename, length, argument_stat)
        if (argument_stat /= 0) then
            call error_exit("create_filename_from_argument: error")
        end if
        call check_file_exists(filename(1:length))
        filename_i = filename(1:length)
    end subroutine create_filename_from_argument

end module procedures_command_arguments
