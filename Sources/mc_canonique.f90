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
    real(DP) :: activExInv, activExInvSum ! inverse de l'activité
    integer, parameter :: nWidom = sph_Ncol ! nombre de particules test
    real(DP) :: tIni, tFin
    
    integer, parameter :: unitObs = 10, unitSnapIni = 11, unitSnapFin = 12, &
        unitRapport = 13, unitObsTherm = 14, unit_dx = 15, unitSnapEnCours = 16
        
    type(Component) :: sph
    
    sph = sph_constructor()
        
    write(*, *) "MC_C+Neigh : Vol =", product(Lsize)
    
    Nrejects = 0
    tauxRejectsSum = 0._DP
    enTotSum = 0._DP
    activExInvSum = 0._DP
    
    open(unit=unitRapport, recl=4096, file="rapport.out", status='new', &
        action='write')         ! contre line folding
    call sph%rapport(nWidom, unitRapport)
    call init_random_seed(unitRapport)
    
    ! Initial condition
    call condIni(unitRapport, sph%X)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call sph%snapShot(unitSnapIni)
    close(unitSnapIni)
    call sph%overlapTest()
    enTot = sph%enTotCalc()
    call sph%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Début des cycles"
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
        
            call sph%adapt_dx(iStep, tauxRejectsSum, unitRapport)
            write(unit_dx, *) iStep, sph%getDx(), tauxRejectsSum/real(iStep, DP)
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
    write(*, *) "Fin des cycles"

! End -----------------------------------------------------

    call sph%overlapTest()
    
    write(unitRapport, *) "Test de consistence :"
    write(unitRapport, *) "    enTot_mc_c = ", enTot
    write(unitRapport, *) "    enTot_calc = ", sph%enTotCalc()
    
    call mcResults(enTotSum, activExInvSum, tauxRejectsSum, tFin-tIni,&
        unitRapport)
    close(unitRapport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call sph%snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call sph%destructor()
    
end program mc_cano
