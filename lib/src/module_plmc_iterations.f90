module module_plmc_iterations

use json_module, only: json_file
use procedures_checks, only: check_data_found, check_positive

implicit none

private
public :: num_tuning_steps, num_steps, plmc_set_num_steps

    integer, protected :: num_tuning_steps, num_steps

contains

    subroutine plmc_set_num_steps(input_data)
        type(json_file), intent(inout) :: input_data

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = "Monte Carlo.number of tuning steps"
        call input_data%get(data_field, num_tuning_steps, data_found)
        call check_data_found(data_field, data_found)
        call check_positive("plmc_set_num_steps", "num_tuning_steps", num_tuning_steps)
        data_field = "Monte Carlo.number of steps"
        call input_data%get(data_field, num_steps, data_found)
        call check_data_found(data_field, data_found)
        call check_positive("plmc_set_num_steps", "num_steps", num_steps)
        deallocate(data_field)
    end subroutine plmc_set_num_steps

end module module_plmc_iterations
