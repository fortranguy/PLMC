module procedures_plmc_help

use, intrinsic :: iso_fortran_env, only: output_unit
use data_strings, only: max_word_length
use procedures_errors, only: error_exit

implicit none

private
public :: plmc_catch_generating_help, plmc_catch_exploring_help, plmc_catch_density_help, &
    plmc_catch_radial_help

contains

    subroutine plmc_catch_generating_help()
        call plmc_catch_help_core("plmc_catch_generating_help", &
            "Please provide generating.json.")
    end subroutine plmc_catch_generating_help

    subroutine plmc_catch_exploring_help()
        call plmc_catch_help_core("plmc_catch_exploring_help", &
            "Please provide generating.json, exploring.json and snap shots.")
    end subroutine plmc_catch_exploring_help

    subroutine plmc_catch_density_help()
        call plmc_catch_help_core("plmc_catch_density_help", &
            "Please provide generating.json, exploring.json and snap shots.")
    end subroutine plmc_catch_density_help

    subroutine plmc_catch_radial_help()
        call plmc_catch_help_core("plmc_catch_radial_help", &
            "Please provide generating.json, exploring.json and snap shots.")
    end subroutine plmc_catch_radial_help

    subroutine plmc_catch_help_core(procedure_name, message)
        character(len=*), intent(in) :: procedure_name, message

        character(len=max_word_length) :: argument
        integer :: length, status

        call get_command_argument(1, argument, length, status)
        if (status /= 0) call error_exit("in procedures_plmc_help: "//procedure_name//&
            ". Type -h for help.")

        select case (argument)
            case ("-h", "--help")
                write(output_unit, *) message
                stop
        end select
    end subroutine plmc_catch_help_core

end module procedures_plmc_help
