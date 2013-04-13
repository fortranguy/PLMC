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
    type(Interacting) :: inter_sph
    type(Observables) :: inter_obs
    
    real(DP) :: tIni, tFin
    
    type(FileUnits) :: units
        
    inter_sph = inter_constructor()
    call inter_obs%init()
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
        action='write')
    call inter_sph%report(unitReport)
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call initialCondition(unitReport, inter_sph%X)
    open(unit=unitSnapIni, recl=4096, file="snapShotIni.out", status='new', &
        action='write')
        call inter_sph%snapShot(unitSnapIni)
    close(unitSnapIni)
    
    call inter_sph%overlapTest()
    inter_obs%ePot_total = inter_sph%ePot_total()
    call inter_sph%cols_to_cells()
    
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
            call inter_sph%mcMove(inter_obs%ePot_total, inter_obs%Nrejects)
        end do
        
        call inter_sph%widom(inter_obs%activExInv)
        
        if (iStep <= Ntherm) then
        
            call inter_sph%adapt_dx(iStep, inter_obs%rejectsRateSum, unitReport)
            write(unit_dx, *) iStep, inter_sph%get_dx(), &
                inter_obs%rejectsRateSum/real(iStep, DP)
            write(unitObsTherm, *) iStep, inter_obs%ePot_total, &
                inter_obs%activExInv
        
        else
            
            call inter_obs%add()
            write(unitObs, *) iStep, inter_obs%ePot_total, inter_obs%activExInv
            
            if (snap) then
                call inter_sph%snapShot(unitSnapShots)
            end if
            
        end if
        
        call inter_obs%addReject()
    
    end do
    call cpu_time(tFin)
    
    close(unitSnapShots)
    close(unit_dx)
    close(unitObsTherm)
    close(unitObs)
    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call inter_sph%overlapTest()
    call inter_sph%consistTest(inter_obs%ePot_total, unitReport)
    
    call inter_obs%results(tFin-tIni, unitReport)
    close(unitReport)
    
    open(unit=unitSnapFin, recl=4096, file="snapShotFin.out", status='new', &
        action='write')
        call inter_sph%snapShot(unitSnapFin)
    close(unitSnapFin)
    
    call inter_sph%destructor()
    
end program mc_canonical
