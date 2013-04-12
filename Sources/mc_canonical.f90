!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_interacting
use class_observables

implicit none

! Beginning

    ! Initialisation
    
    integer :: iStep, iMove
    type(Interacting) :: inter
    type(Observables) :: interObs
    
    real(DP) :: tIni, tFin
    
    integer, parameter :: unitObs = 10, unitSnapIni = 11, unitSnapFin = 12, &
        unitReport = 13, unitObsTherm = 14, unit_dx = 15, unitSnapShots = 16
        
    inter = inter_constructor()
    call interObs%init()
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
        action='write')
    call inter%report(unitReport)
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call initialCondition(unitReport, inter%X)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call inter%snapShot(unitSnapIni)
    close(unitSnapIni)
    
    call inter%overlapTest()
    interObs%ePot_total = inter%ePot_total()
    call inter%cols_to_cells()
    
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
            call inter%mcMove(interObs%ePot_total, interObs%Nrejects)
        end do
        
        call inter%widom(interObs%activExInv)
        
        if (iStep <= Ntherm) then
        
            call inter%adapt_dx(iStep, interObs%rejectsRateSum, unitReport)
            write(unit_dx, *) iStep, inter%get_dx(), &
                interObs%rejectsRateSum/real(iStep, DP)
            write(unitObsTherm, *) iStep, interObs%ePot_total, &
                interObs%activExInv
        
        else
            
            call interObs%add()
            write(unitObs, *) iStep, interObs%ePot_total, interObs%activExInv
            
            if (snap) then
                call inter%snapShot(unitSnapShots)
            end if
            
        end if
        
        call interObs%addRej()
    
    end do
    call cpu_time(tFin)
    
    close(unitSnapShots)
    close(unit_dx)
    close(unitObsTherm)
    close(unitObs)
    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call inter%overlapTest()
    call inter%consistTest(interObs%ePot_total, unitReport)
    
    call mcResults(interObs, tFin-tIni, unitReport)
    close(unitReport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call inter%snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call inter%destructor()
    
end program mc_canonical
