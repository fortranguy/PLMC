module procedures_plmc_help

use, intrinsic :: iso_fortran_env, only: output_unit
use data_strings, only: max_word_length
use procedures_errors, only: error_exit

implicit none

private
public :: plmc_catch_help, plmc_widom_catch_help

contains

    subroutine plmc_catch_help()
        call plmc_catch_help_core("Please provide a data file in json format.")
    end subroutine plmc_catch_help

    subroutine plmc_widom_catch_help()
        call plmc_catch_help_core("Please provide snaps, input_data.json and post.json.")
    end subroutine plmc_widom_catch_help

    subroutine plmc_catch_help_core(message)
        character(len=*), intent(in) :: message

        character(len=max_word_length) :: argument
        integer :: length, status

        call get_command_argument(1, argument, length, status)
        if (status /= 0) call error_exit("in procedures_plmc_help: plmc_help(). Type -h for help.")

        select case (argument)
            case ("-h", "--help")
                write(output_unit, *) message
                stop
        end select
    end subroutine plmc_catch_help_core

end module procedures_plmc_help
