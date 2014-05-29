!> \brief Bunching program to assess the error of an observable

program bunching

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use data_precisions, only: DP
use data_monte_carlo, only: num_equilibrium_steps

implicit none

    ! Physics
    
    integer :: Nobs
    integer :: iBunching, Nbunching
    integer :: i_step, iStepIn, NstepVar
    real(DP), dimension(:), allocatable :: sumVal, sumValSqr
    real(DP), dimension(:), allocatable :: error

    ! Numerical
    
    integer :: obs_unit, bunch_unit
    character(len=1) :: comment_symbol
    integer, parameter :: NdataMi = num_equilibrium_steps/2
    real(DP), dimension(:, :), allocatable :: dataIn
    real(DP), dimension(:, :), allocatable :: dataOut
    
    character(len=4096) :: file
    integer :: length, stat
    logical :: exist
    
    call get_command_argument(1, file, length, stat)
    if (stat /= 0) error stop "error get_command_argument"
    inquire(file=file(1:length), exist=exist)
    if (.not.exist) then
        write(error_unit, *) "missing file: ", file(1:length)
        error stop
    end if
    
    open(newunit=obs_unit, recl=4096, file=file(1:length), status='old', action='read')
    read(obs_unit, *) comment_symbol, Nobs
    write(output_unit, *) "Nobs = ", Nobs
    
    allocate(sumVal(Nobs))
    allocate(sumValSqr(Nobs))
    allocate(error(Nobs))
    
    allocate(dataIn(Nobs, 2*NdataMi))
    allocate(dataOut(Nobs, NdataMi))
    
    write(output_unit, *) "num_equilibrium_steps = ", num_equilibrium_steps
    Nbunching = int(log(real(num_equilibrium_steps, DP))/log(2._DP))
    write(output_unit, *) "Nbunching = ", Nbunching

    NstepVar = num_equilibrium_steps
    
    open(newunit=bunch_unit, recl=4096, file=file(1:length-4)//"_bunched.out", action='write')
    
    do iBunching = 1, Nbunching
    
        NstepVar = NstepVar/2
        
        ! Read
        if (iBunching == 1) then
            do i_step = 1, 2*NstepVar
                read(obs_unit, *) iStepIn, dataIn(:, i_step)
            end do
        else
            do i_step = 1, 2*NstepVar
                dataIn(:, i_step) = dataOut(:, i_step)
            end do
        end if
    
        ! Bunch
        sumVal = 0.
        sumValSqr = 0.
        
        do i_step = 1, NstepVar
            sumVal(:) = sumVal(:) + dataIn(:, 2*i_step-1) + dataIn(:, 2*i_step)
            sumValSqr(:) = sumValSqr(:) + dataIn(:, 2*i_step-1)**2 + dataIn(:, 2*i_step)**2
            dataOut(:, i_step) = (dataIn(:, 2*i_step-1) + dataIn(:, 2*i_step))/2.
        end do
        
        error(:) = sumValSqr(:)/real(2*NstepVar, DP) - (sumVal(:)/real(2*NstepVar, DP))**2
        error(:) = sqrt(error(:)/real(2*NstepVar, DP))
        
        ! Results
        write(bunch_unit, *) iBunching, error(:)
        
    end do
    
    close(bunch_unit)
    
    deallocate(dataOut)
    deallocate(dataIn)
    
    deallocate(error)
    deallocate(sumValSqr)
    deallocate(sumVal)
    
    close(obs_unit)

end program bunching
