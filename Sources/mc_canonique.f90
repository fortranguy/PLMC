program mc_cano

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_component

implicit none

! Beginning ---------------------------------------------------

    ! Initialisation
    
    integer :: iStep, iMove
    integer :: Nrejects
    real(DP) :: tauxRejectsSum
    real(DP) :: enTot, enTotSum
    !	inverse of activity :
    real(DP) :: activExInv, activExInvSum
    !	number of test particles :
    integer, parameter :: nWidom = sph_Ncol
    real(DP) :: tIni, tFin
    type(Component) :: sph
    
    integer, parameter :: unitObs = 10, unitSnapIni = 11, unitSnapFin = 12, &
        unitReport = 13, unitObsTherm = 14, unit_dx = 15, unitSnapEnCours = 16
        
    sph = sph_constructor()
    
    Nrejects = 0
    tauxRejectsSum = 0._DP
    enTotSum = 0._DP
    activExInvSum = 0._DP
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
        action='write')
    call sph%report(nWidom, unitReport)
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call condIni(unitReport, sph%X)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call sph%snapShot(unitSnapIni)
    close(unitSnapIni)
    
    call sph%overlapTest()
    enTot = sph%enTotCalc()
    call sph%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Beginning of cycles"
    open(unit=unitObs, recl=4096, file="obs.out", status='new', &
        action='write')
    open(unit=unitObsTherm, recl=4096, file="obsTherm.out", status='new', &
        action='write')
    open(unit=unit_dx, recl=4096, file="dx.out", status='new', &
        action='write')
    open(unit=unitSnapEnCours, recl=4096, file="snap.shot", status='new', &
        action='write')
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove   
            call sph%mcMove(enTot, Nrejects)
        end do
        
        call sph%widom(nWidom, activExInv)
        
        if (iStep <= Ntherm) then
        
            call sph%adapt_dx(iStep, tauxRejectsSum, unitReport)
            write(unit_dx, *) iStep, sph%get_dx(), tauxRejectsSum/real(iStep, DP)
            write(unitObsTherm, *) iStep, enTot, activExInv
        
        else
            
            enTotSum = enTotSum + enTot
            activExInvSum = activExInvSum + activExInv
            write(unitObs, *) iStep, enTot, activExInv
            
            if (snap) then
                call sph%snapShot(unitSnapEnCours)
            end if
            
        end if
        
        tauxRejectsSum = tauxRejectsSum + real(Nrejects, DP)/real(Nmove, DP)
        Nrejects = 0
    
    end do
    call cpu_time(tFin)
    
    close(unitSnapEnCours)
    close(unit_dx)
    close(unitObsTherm)
    close(unitObs)
    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call sph%overlapTest()
    
    write(unitReport, *) "Consistency test:"
    write(unitReport, *) "    enTot_mc_c = ", enTot
    write(unitReport, *) "    enTot_calc = ", sph%enTotCalc()
    
    call mcResults(enTotSum, activExInvSum, tauxRejectsSum, tFin-tIni,&
        unitReport)
    close(unitReport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call sph%snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call sph%destructor()
    
end program mc_cano
