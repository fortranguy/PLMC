!> \brief Bunching program to assess the error of an observable

program bunching

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_monteCarlo, only : Nstep

implicit none

    ! Physics 
    
    integer, parameter :: nObs = 2
    integer :: iBunching, nBunching
    integer :: iStep, iStepIn, NstepVar
    real(DP), dimension(nObs) :: sumVal, sumValSqr
    real(DP), dimension(nObs) :: error

    ! Numerical    
    
    integer :: obs_unit
    integer, parameter :: nDataMi = Nstep/2
    real(DP), dimension(nObs, 2*nDataMi) :: dataIn
    real(DP), dimension(nObs, nDataMi) :: dataOut
    
    character(len=4096) :: file_name
    integer :: length, file_stat
    
    call get_command_argument(1, file_name, length, file_stat)
    if (file_stat /= 0) stop "error get_command_argument"
    
    nBunching = int(log(real(Nstep, DP))/log(2._DP))
    write(output_unit, *) "nBunching = ", nBunching

    NstepVar = Nstep
    
    write(output_unit, *)
    
    open(unit=11, file="bunching_eTot.out")
    open(unit=12, file="bunching_activInv.out")
    
    do iBunching = 1, nBunching
    
        write(output_unit, *) "iBunching = ", iBunching, "NstepVar = ", NstepVar
        NstepVar = NstepVar/2
        
        ! Read
        
        if (iBunching == 1) then
        
            open(newunit=obs_unit, recl=4096, file=file_name(1:length), status='old', action='read')
            do iStep = 1, 2*NstepVar
                read(obs_unit, *) iStepIn, dataIn(1, iStep), dataIn(2, iStep)
            end do
            close(obs_unit)
            
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
        
        write(11, *) iBunching, sumVal(1)/real(2*NstepVar, DP), error(1)
        write(12, *) iBunching, sumVal(2)/real(2*NstepVar, DP), error(2)
        
    end do

    close(12)
    close(11)

end program bunching
