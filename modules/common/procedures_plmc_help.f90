module procedures_plmc_help

use, intrinsic :: iso_fortran_env, only: output_unit
use data_strings, only: max_word_length
use procedures_errors, only: error_exit

implicit none

private
public :: plmc_catch_help

contains

    subroutine plmc_catch_help()

        character(len=max_word_length) :: argument
        integer :: length, status

        call get_command_argument(1, argument, length, status)
        if (status /= 0) call error_exit("in procedures_plmc_help: plmc_help(). Type -h for help.")

        select case (argument)
            case ("-h", "--help")
                write(output_unit, *) "Please provide a data file in json format."
                stop
        end select
    end subroutine plmc_catch_help

end module procedures_plmc_help
