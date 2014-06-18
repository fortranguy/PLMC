!> \brief Histogram program of observables: energy

program histogram

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use json_module, only: json_file, json_initialize
use module_data, only: test_data_found
use module_post_processing_arguments, only: argument_to_file

implicit none

    integer :: num_observables
    integer :: num_steps, i_step, i_step_in
    integer :: i_obsersable, i_distribution
    real(DP), dimension(:, :), allocatable :: limit_values
    integer, dimension(:, :), allocatable :: num_distribution
    real(DP), dimension(:), allocatable :: energy_distribution_function
    real(DP) :: energy_i, energy_delta
    
    integer :: observables_unit, histogram_unit
    character(len=1) :: comment_symbol
    real(DP), dimension(:, :), allocatable :: observables
    
    type(json_file) :: json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: file_name
    integer :: length
    
    call json_initialize()
    call json%load_file(filename = "data.json")
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call json%get(data_name, num_steps, found)
    call test_data_found(data_name, found)
    
    data_name = "Histogram.energy delta"
    call json%get(data_name, energy_delta, found)
    call test_data_found(data_name, found)
    
    call json%destroy()
    
    call argument_to_file(1, file_name, length)  
    
    open(newunit=observables_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    read(observables_unit, *) comment_symbol, num_observables
    write(output_unit, *) "num_observables = ", num_observables
    
    allocate(observables(num_observables, num_steps))
    allocate(limit_values(2, num_observables))
    allocate(num_distribution(2, num_observables))
    
    write(output_unit, *) "num_steps = ", num_steps
    
    do i_step = 1, num_steps
        read(observables_unit, *) i_step_in, observables(:, i_step)
    end do
    
    do i_obsersable = 1, num_observables
        limit_values(1, i_obsersable) = minval(observables(i_obsersable, :))
        limit_values(2, i_obsersable) = maxval(observables(i_obsersable, :))
    end do
    num_distribution(:, 1) = int(limit_values(:, 1)/energy_delta)
    allocate(energy_distribution_function(num_distribution(1, 1): num_distribution(2, 1)))
    
    energy_distribution_function(:) = 0._DP
    do i_step = 1, num_steps
        i_distribution = int(observables(1, i_step)/energy_delta)
        energy_distribution_function(i_distribution) = &
            energy_distribution_function(i_distribution) + 1._DP
    end do
    energy_distribution_function(:) = energy_distribution_function(:) /  real(num_steps, DP) / &
                                      energy_delta
    
    open(newunit=histogram_unit, recl=4096, file=file_name(1:length-4)//"_energy_histogram.out", &
         action='write')
    do i_distribution = num_distribution(1, 1), num_distribution(2, 1)
        energy_i = (real(i_distribution, DP)+0.5_DP) * energy_delta
        write(histogram_unit, *) energy_i, energy_distribution_function(i_distribution)
    end do
    
    close(histogram_unit)
    close(observables_unit)
    
    deallocate(energy_distribution_function)
    deallocate(num_distribution)
    deallocate(limit_values)
    deallocate(observables)

end program histogram
