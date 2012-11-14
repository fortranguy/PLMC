program bunching

use data_constants
use data_mc

implicit none

    ! Physique 
    
    integer, parameter :: nObs = 2
    character(len=20) :: nBunchingStr
    integer :: iBunching, nBunching
    integer :: iStep, iStepIn, NstepVar
    real(DP), dimension(nObs) :: sumVal, sumValSqr
    real(DP), dimension(nObs) :: error

    ! Numérique    
    
    integer :: longueur, statut
    integer, parameter :: nDataMi = 2**22 ! Attention
    real(DP), dimension(nObs, 2*nDataMi) :: dataIn
    real(DP), dimension(nObs, nDataMi) :: dataOut
    
    call get_command_argument(1, nBunchingStr, longueur, statut)
    if (statut /= 0) stop "erreur get_command_argument"
    if (command_argument_count() > 1) stop "Trop d'arguments"
    read(nBunchingStr, *) nBunching
    write(*, *) "nBunching = ", nBunching

    NstepVar = Nstep
    
    write(*,*)
    
    open(unit=11, file="bunching_eTot.out")
    open(unit=12, file="bunching_activInv.out")
    
    do iBunching = 1, nBunching
    
        write(*,*) "iBunching = ", iBunching, "NstepVar = ", NstepVar    
        NstepVar = NstepVar/2
        
        ! Lecture
        
        if (iBunching == 1) then
        
            open(unit=10, file="obs.out")
            do iStep = 1, 2*NstepVar
                read(10, *) iStepIn, dataIn(1, iStep), dataIn(2, iStep)
            end do
            close(10)
            
        else
        
            do iStep = 1, 2*NstepVar
                dataIn(:, iStep) = dataOut(:, iStep)
            end do
            
        end if
    
        ! Bunching
        
        sumVal = 0.
        sumValSqr = 0.
        
        do iStep = 1, NstepVar
        
            sumVal(:) = sumVal(:) + dataIn(:, 2*iStep-1) + dataIn(:, 2*iStep)
            sumValSqr(:) = sumValSqr(:) + dataIn(:, 2*iStep-1)**2 + &
                dataIn(:, 2*iStep)**2
            dataOut(:, iStep) = (dataIn(:, 2*iStep-1) + dataIn(:, 2*iStep))/2.
            
        end do
        
        error(:) = sumValSqr(:)/real(2*NstepVar, DP) - &
            (sumVal(:)/real(2*NstepVar, DP))**2
        error(:) = sqrt(error(:)/real(2*NstepVar, DP))
        
        ! Résultats
        
        write(*, *) "moyennes = ", sumVal(:)/real(2*NstepVar, DP)
        write(*, *) "erreurs = ", error(:)
        write(11, *) iBunching, sumVal(1)/real(2*NstepVar, DP), error(1)
        write(12, *) iBunching, sumVal(2)/real(2*NstepVar, DP), error(2)
        
    end do

    close(12)
    close(11)

end program bunching
