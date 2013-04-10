!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_component

implicit none

! Beginning

    ! Initialisation
    
    integer :: iStep, iMove
    integer :: Nrejects
    real(DP) :: tauxRejectsSum
    real(DP) :: enTot, enTotSum
    real(DP) :: activExInv, activExInvSum !< inverse of activity
    integer, parameter :: nWidom = inter_Ncol !< number of test particles
    real(DP) :: tIni, tFin
    type(Interacting) :: inter
    
    integer, parameter :: unitObs = 10, unitSnapIni = 11, unitSnapFin = 12, &
        unitReport = 13, unitObsTherm = 14, unit_dx = 15, unitSnapShots = 16
        
    inter = inter_constructor()
    
    Nrejects = 0
    tauxRejectsSum = 0._DP
    enTotSum = 0._DP
    activExInvSum = 0._DP
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
        action='write')
    call intr%report(nWidom, unitReport)
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call initialCondition(unitReport, intr%X)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call intr%snapShot(unitSnapIni)
    close(unitSnapIni)
    
    call intr%overlapTest()
    enTot = intr%enTotCalc()
    call intr%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Beginning of cycles"
    open(unit=unitObs, recl=4096, file="obs.out", status='new', &
        action='write')
    open(unit=unitObsTherm, recl=4096, file="obsTherm.out", status='new', &
        action='write')
    open(unit=unit_dx, recl=4096, file="dx.out", status='new', &
        action='write')
    open(unit=unitSnapShots, recl=4096, file="snap.shot", status='new', &
        action='write')
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove   
            call intr%mcMove(enTot, Nrejects)
        end do
        
        call intr%widom(nWidom, activExInv)
        
        if (iStep <= Ntherm) then
        
            call intr%adapt_dx(iStep, tauxRejectsSum, unitReport)
            write(unit_dx, *) iStep, intr%get_dx(), tauxRejectsSum/real(iStep, DP)
            write(unitObsTherm, *) iStep, enTot, activExInv
        
        else
            
            enTotSum = enTotSum + enTot
            activExInvSum = activExInvSum + activExInv
            write(unitObs, *) iStep, enTot, activExInv
            
            if (snap) then
                call intr%snapShot(unitSnapShots)
            end if
            
        end if
        
        tauxRejectsSum = tauxRejectsSum + real(Nrejects, DP)/real(Nmove, DP)
        Nrejects = 0
    
    end do
    call cpu_time(tFin)
    
    close(unitSnapShots)
    close(unit_dx)
    close(unitObsTherm)
    close(unitObs)
    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call intr%overlapTest()
    
    write(unitReport, *) "Consistency test:"
    write(unitReport, *) "    enTot_mc_c = ", enTot
    write(unitReport, *) "    enTot_calc = ", intr%enTotCalc()
    
    call mcResults(enTotSum, activExInvSum, tauxRejectsSum, tFin-tIni,&
        unitReport)
    close(unitReport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call intr%snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call intr%destructor()
    
end program mc_canonical
