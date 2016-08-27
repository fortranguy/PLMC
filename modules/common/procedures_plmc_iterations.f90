module procedures_plmc_iterations

use json_module, only: json_file
use data_input_prefixes, only: changes_prefix, writers_prefix
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found, check_positive
use procedures_writers_inquirers, only: write_coordinates

implicit none

private
public :: plmc_set_num_steps, plmc_set_num_snaps

contains

    subroutine plmc_set_num_steps(num_tuning_steps, num_steps, generating_data)
        integer, intent(out) :: num_tuning_steps, num_steps
        type(json_file), intent(inout) :: generating_data

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = changes_prefix//"number of tuning steps"
        call generating_data%get(data_field, num_tuning_steps, data_found)
        call check_data_found(data_field, data_found)
        call check_positive("plmc_set_num_steps", "num_tuning_steps", num_tuning_steps)
        data_field = "Monte Carlo.number of steps"
        call generating_data%get(data_field, num_steps, data_found)
        call check_data_found(data_field, data_found)
        call check_positive("plmc_set_num_steps", "num_steps", num_steps)
    end subroutine plmc_set_num_steps

    subroutine plmc_set_num_snaps(num_snaps, generating_data)
        integer, intent(out) :: num_snaps
        type(json_file), intent(inout) :: generating_data

        if (.not.write_coordinates(generating_data, writers_prefix)) then
            call error_exit("Coordinates weren't written.")
        end if
        num_snaps = command_argument_count() - 2
        if (num_snaps == 0) call error_exit("No snaps given.")
    end subroutine plmc_set_num_snaps

end module procedures_plmc_iterations
