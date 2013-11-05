!> \brief Bunching program to assess the error of an observable

program bunching

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_monteCarlo, only : Nstep

implicit none

    ! Physics 
    
    integer :: Nobs
    integer :: iBunching, Nbunching
    integer :: iStep, iStepIn, NstepVar
    real(DP), dimension(:), allocatable :: sumVal, sumValSqr
    real(DP), dimension(:), allocatable :: error

    ! Numerical    
    
    integer :: obs_unit, bunch_unit
    character(len=1) :: comment_symbol
    integer, parameter :: NdataMi = Nstep/2
    real(DP), dimension(:, :), allocatable :: dataIn
    real(DP), dimension(:, :), allocatable :: dataOut
    
    character(len=4096) :: file_name
    integer :: length, file_stat
    
    call get_command_argument(1, file_name, length, file_stat)
    if (file_stat /= 0) stop "error get_command_argument"
    
    open(newunit=obs_unit, recl=4096, file=file_name(1:length), status='old', action='read')
    read(obs_unit, *) comment_symbol, Nobs
    
    allocate(sumVal(Nobs))
    allocate(sumValSqr(Nobs))
    allocate(error(Nobs))
    
    allocate(dataIn(Nobs, 2*NdataMi))
    allocate(dataOut(Nobs, NdataMi))
    
    Nbunching = int(log(real(Nstep, DP))/log(2._DP))
    write(output_unit, *) "Nbunching = ", Nbunching

    NstepVar = Nstep
    
    open(newunit=bunch_unit, recl=4096, file=file_name(1:length-4)//"_bunched.out", action='write')
    
    do iBunching = 1, Nbunching
    
        write(output_unit, *) "iBunching = ", iBunching, "NstepVar = ", NstepVar
        NstepVar = NstepVar/2
        
        ! Read
        
        if (iBunching == 1) then        
            
            do iStep = 1, 2*NstepVar
                read(obs_unit, *) iStepIn, dataIn(:, iStep)
            end do            
            
        else
        
            do iStep = 1, 2*NstepVar
                dataIn(:, iStep) = dataOut(:, iStep)
            end do
            
        end if
    
        ! Bunch
        
        sumVal = 0.
        sumValSqr = 0.
        
        do iStep = 1, NstepVar
        
            sumVal(:) = sumVal(:) + dataIn(:, 2*iStep-1) + dataIn(:, 2*iStep)
            sumValSqr(:) = sumValSqr(:) + dataIn(:, 2*iStep-1)**2 + dataIn(:, 2*iStep)**2
            dataOut(:, iStep) = (dataIn(:, 2*iStep-1) + dataIn(:, 2*iStep))/2.
            
        end do
        
        error(:) = sumValSqr(:)/real(2*NstepVar, DP) - (sumVal(:)/real(2*NstepVar, DP))**2
        error(:) = sqrt(error(:)/real(2*NstepVar, DP))
        
        ! Results
        
        write(bunch_unit, *) iBunching, error(:)
        
    end do
    
    deallocate(sumVal)
    deallocate(sumValSqr)
    deallocate(error)
    
    deallocate(dataIn)
    deallocate(dataOut)
    
    close(obs_unit)
    close(bunch_unit)

end program bunching
