!> \brief Bunching program to assess the error of an observable

program bunching

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: data_filename, &
                       test_file_exists, test_data_found
use module_arguments, only: arg_to_file

implicit none

    ! Physics
    
    integer :: num_steps
    integer :: i_step, i_step_in
    integer :: num_observables
    integer :: i_bunching, num_bunching
    real(DP), dimension(:), allocatable :: values_sum, values_sqr_sum
    real(DP), dimension(:), allocatable :: error

    ! Numerical
    
    integer :: observables_unit, bunching_unit
    character(len=1) :: comment_symbol
    real(DP), dimension(:, :), allocatable :: data_in
    real(DP), dimension(:, :), allocatable :: data_out
    
    type(json_file) :: data_json
    character(len=4096) :: data_name
    logical :: found
    character(len=4096) :: filename
    integer :: length
    
    call json_initialize()
    
    call test_file_exists(data_filename)
    call data_json%load_file(filename = data_filename)
    
    data_name = "Monte Carlo.number of equilibrium steps"
    call data_json%get(data_name, num_steps, found)
    call test_data_found(data_name, found)
    
    call data_json%destroy()
    
    call arg_to_file(1, filename, length)
    open(newunit=observables_unit, recl=4096, file=filename(1:length), status='old', action='read')
    read(observables_unit, *) comment_symbol, num_observables
    write(output_unit, *) "num_observables = ", num_observables
    
    allocate(values_sum(num_observables))
    allocate(values_sqr_sum(num_observables))
    allocate(error(num_observables))
    
    allocate(data_in(num_observables, num_steps))
    allocate(data_out(num_observables, num_steps/2))
    
    write(output_unit, *) "num_steps = ", num_steps
    num_bunching = int(log(real(num_steps, DP))/log(2._DP))
    write(output_unit, *) "num_bunching = ", num_bunching
    
    open(newunit=bunching_unit, recl=4096, file=filename(1:length-4)//"_bunching.out", action='write')
    
    do i_bunching = 1, num_bunching
    
        num_steps = num_steps/2
        
        ! Read
        if (i_bunching == 1) then
            do i_step = 1, 2*num_steps
                read(observables_unit, *) i_step_in, data_in(:, i_step)
            end do
        else
            do i_step = 1, 2*num_steps
                data_in(:, i_step) = data_out(:, i_step)
            end do
        end if
    
        ! Bunch
        values_sum = 0.
        values_sqr_sum = 0.
        
        do i_step = 1, num_steps
            values_sum(:) = values_sum(:) + data_in(:, 2*i_step-1) + data_in(:, 2*i_step)
            values_sqr_sum(:) = values_sqr_sum(:) + data_in(:, 2*i_step-1)**2 + data_in(:, 2*i_step)**2
            data_out(:, i_step) = (data_in(:, 2*i_step-1) + data_in(:, 2*i_step))/2.
        end do
        
        error(:) = values_sqr_sum(:)/real(2*num_steps, DP) - (values_sum(:)/real(2*num_steps, DP))**2
        error(:) = sqrt(error(:)/real(2*num_steps, DP))
        
        ! Results
        write(bunching_unit, *) i_bunching, error(:)
        
    end do
    
    close(bunching_unit)
    
    deallocate(data_out)
    deallocate(data_in)
    
    deallocate(error)
    deallocate(values_sqr_sum)
    deallocate(values_sum)
    
    close(observables_unit)

end program bunching
