!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_hardSpheres
use class_interactingSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none

! Beginning

    ! Initialisation
    
    integer :: iStep, iMove
    integer :: iColRand
    real(DP) :: rand
    real(DP) :: tIni, tFin
    
    integer, parameter :: unitReport = 10
    
    ! Hard spheres
    type(HardSpheres) :: hard_sph
    type(Observables) :: hard_obs
    type(Units) :: hard_io
    
    ! Interacting spheres
    type(InteractingSpheres) :: inter_sph !< Monte-Carlo subroutines
    type(Observables) :: inter_obs !< e.g. Energy
    type(Units) :: inter_io        !< input/output files
    
    call hard_sph%construct()
    call hard_obs%init()
    call hard_io%open("hard")
        
    call inter_sph%construct()
    call inter_obs%init()
    call inter_io%open("inter")
    
    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)    
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
    	action='write')
	call report(unitReport)
    call hard_sph%report(hard_io%report)
    call inter_sph%report(inter_io%report)
    
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call initialCondition(hard_sph%X, hard_io%report)
    call hard_sph%overlapTest()
    hard_obs%ePot_total = 0._DP
    call hard_sph%snapShot(hard_io%snapIni)
    call hard_sph%cols_to_cells()
    
    call initialCondition(inter_sph%X, inter_io%report)
    call inter_sph%overlapTest()
    inter_obs%ePot_total = inter_sph%ePot_total()
    call inter_sph%snapShot(inter_io%snapIni)
    call inter_sph%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove
        
            call random_number(rand)
            iColRand = int(rand*real(hard_Ncol+inter_Ncol, DP)) + 1            
            if (iColRand<=hard_Ncol) then
                call hard_sph%move(hard_obs%Nrej)
                hard_obs%Nmove = hard_obs%Nmove + 1
            else
                call inter_sph%move(inter_obs%ePot_total, inter_obs%Nrej)
                inter_obs%Nmove = inter_obs%Nmove + 1
            end if            
            
        end do
        
        call hard_sph%widom(hard_obs%activExInv)
        call inter_sph%widom(inter_obs%activExInv)
        
        if (iStep <= Ntherm) then
        
            call hard_sph%adaptDx(iStep, hard_obs%rejRateSum, &
                hard_io%report)
            write(hard_io%dx, *) iStep, hard_sph%getDx(), &
                hard_obs%rejRateSum/real(iStep, DP)
            write(hard_io%obsTherm, *) iStep, hard_obs%ePot_total, &
                hard_obs%activExInv
        
            call inter_sph%adaptDx(iStep, inter_obs%rejRateSum, &
                inter_io%report)
            write(inter_io%dx, *) iStep, inter_sph%getDx(), &
                inter_obs%rejRateSum/real(iStep, DP)
            write(inter_io%obsTherm, *) iStep, inter_obs%ePot_total, &
                inter_obs%activExInv
        
        else
        
            call hard_obs%addPhysical()
            write(hard_io%obs, *) iStep, hard_obs%ePot_total, &
                hard_obs%activExInv
            
            call inter_obs%addPhysical()
            write(inter_io%obs, *) iStep, inter_obs%ePot_total, &
                inter_obs%activExInv
            
            if (snap) then
                call hard_sph%snapShot(hard_io%snapShots)
                call inter_sph%snapShot(inter_io%snapShots)
            end if
            
        end if
        
        call hard_obs%addReject()
        call inter_obs%addReject()
    
    end do
    call cpu_time(tFin)

    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call hard_sph%overlapTest()
    call hard_sph%snapShot(hard_io%snapFin)
    call hard_obs%results(tFin-tIni, hard_io%report)

    call inter_sph%overlapTest()
    call inter_sph%consistTest(inter_obs%ePot_total, inter_io%report)
    call inter_sph%snapShot(inter_io%snapFin)
    call inter_obs%results(tFin-tIni, inter_io%report)
    
    close(unitReport)
    
    call hard_sph%destroy()
    call hard_io%close()
    
    call inter_sph%destroy()
    call inter_io%close()
    
end program mc_canonical
