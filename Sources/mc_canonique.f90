program mc_cano

use data_cell
use data_particles    
use data_mc
use mod_tools
use mod_physique
use data_constants

implicit none

! Début ---------------------------------------------------

    ! Initialisation
    integer :: iStep, iMove
    integer :: Nrejects
    real(DP) :: tauxRejects
    real(DP) :: enTot, enTotSum
    real(DP) :: activExInv, activExInvSum ! inverse de l'activité
    integer, parameter :: nWidom = 800 ! nombre de particules test
    real(DP), parameter :: Lratio = 1._DP ! pour le volume du tirage cf. widom
    
    integer, parameter :: unitObs = 10, unitSnapIni = 11, unitSnapFin = 12, &
        unitRapport = 13, unitObsTherm = 14
    
    Nrejects = 0
    tauxRejects = 0._DP
    enTotSum = 0._DP
    activExInvSum = 0._DP
    call ePotIni()
    !call init_random_seed()
    
    open(unit=unitRapport, recl=4096, file="rapport.out", status='new', &
        action='write')         ! contre line folding
    call rapport(nWidom, Lratio, unitRapport)
    
    ! Condition initiale
    call condIni(unitRapport)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call snapShot(unitSnapIni)
    close(unitSnapIni)
    call overlapTest()
    enTot = enTotCalc()
    
    ! Table des voisins
    call check_CellsSize()
    call alloc_Cells()
    call all_col_to_cell()
    call ini_cell_neighs()
    
! Milieu --------------------------------------------------

    write(*, *) "Début des cycles"
    open(unit=unitObs, recl=4096, file="obs.out", status='new', &
        action='write')
    open(unit=unitObsTherm, recl=4096, file="obsTherm.out", status='new', &
        action='write')
    
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove        
            call mcMove(enTot, Nrejects)    
        end do
        
        call widom(nWidom, Lratio, activExInv)
        
        if (iStep > Ntherm) then
        
            enTotSum = enTotSum + enTot                 
            activExInvSum = activExInvSum + activExInv            
            write(unitObs, *) iStep, enTot, activExInv
            
        else
            
            write(unitObsTherm, *) iStep, enTot, activExInv
            
        end if
        
        tauxRejects = tauxRejects + real(Nrejects, DP)/real(Nmove, DP)
        Nrejects = 0
    
    end do
    
    close(unitObsTherm)
    close(unitObs)
    write(*, *) "Fin des cycles"

! Fin -----------------------------------------------------

    call overlapTest()
    
    call consisTest(enTot, unitRapport)
    call mcResults(enTotSum, activExInvSum, tauxRejects,&
        unitRapport)
    close(unitRapport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call dealloc_Cells()
    
end program mc_cano
