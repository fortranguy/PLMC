!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_mixingPotential
use class_interactingSpheres
use class_hardSpheres
use class_observables
use class_units

implicit none

! Beginning
    
    integer :: iStep, iMove
    integer :: iColRand
    real(DP) :: rand
    real(DP) :: tIni, tFin
    
    real(DP) :: ePot, ePotSum, ePot_total    
    integer, parameter :: unitReport = 10
    integer, parameter :: unitObsTherm = 11, unitObs = 12
    
    ! Type 1 : Interacting spheres
    type(InteractingSpheres) :: type1_sph !< Monte-Carlo subroutines
    type(Observables) :: type1_obs !< e.g. Energy
    type(Units) :: type1_io        !< input/output files
    
    ! Type 2 : Hard spheres
    type(HardSpheres) :: type2_sph
    type(Observables) :: type2_obs
    type(Units) :: type2_io
    
    ! Mixing between 2 types
    type(MixingPotential) :: mix
    real(DP) :: mix_ePot, mix_ePotSum, mix_ePot_total
    integer, parameter :: mix_unitReport = 13
    integer, parameter :: mix_unitObsTherm = 14, mix_unitObs = 15
    
    write(*, *) "Monte-Carlo Mix - Canonical : Volume =", product(Lsize)

    ! Initialisation
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
        action='write')
    open(unit=unitObsTherm, recl=4096, file="obsTherm.out", status='new', &
        action='write')
    open(unit=unitObs, recl=4096, file="obs.out", status='new', &
        action='write')
    
    call type1_sph%construct()
    call type1_obs%init()
    call type1_io%open(type1_sph%getName())
    
    call type2_sph%construct()
    call type2_obs%init()
    call type2_io%open(type2_sph%getName())
    
    call mix%construct()
    mix_ePotSum = 0._DP
    open(unit=mix_unitReport, recl=4096, file="mix_report.out", &
        status='new', action='write')
    open(unit=mix_unitObsTherm, recl=4096, file="mix_obsTherm.out", &
        status='new', action='write')
    open(unit=mix_unitObs, recl=4096, file="mix_obs.out", status='new', &
        action='write')
        
    call report(unitReport)
    call init_random_seed(unitReport)
    
    call type1_sph%report(type1_io%report)
    call type1_sph%printInfo(type1_io%report)
    
    call type2_sph%report(type2_io%report)
    call type2_sph%printInfo(type2_io%report)
    
    call mix%report(mix_unitReport)
    
    ! Initial condition
    
    call initialCondition(type1_sph, type2_sph, mix%getRmin(), unitReport)
    
    call type1_sph%overlapTest()
    type1_obs%ePot = type1_sph%ePot_total()
    call type1_sph%snapShot(type1_io%snapIni)
    call type1_sph%cols_to_cells(type2_sph%X)
    
    call type2_sph%overlapTest()
    type2_obs%ePot = type2_sph%ePot_total()
    call type2_sph%snapShot(type2_io%snapIni)
    call type2_sph%cols_to_cells(type1_sph%X)
    
    call mix%overlapTest(type1_sph%X, type2_sph%X)
    mix_ePot = mix%ePot_total(type1_sph%X, type2_sph%X)
    
    ePot_total = type1_obs%ePot + type2_obs%ePot + mix_ePot
    
! Middle --------------------------------------------------
        
    write(*, *) "Beginning of cycles"    
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove
        
            call random_number(rand)
            iColRand = int(rand*real(Ncol, DP)) + 1            
            if (iColRand <= type1_sph%getNcol()) then
                call type1_sph%move(type2_sph, mix, type1_obs%ePot, &
                    mix_ePot, type1_obs%Nrej)                    
                type1_obs%Nmove = type1_obs%Nmove + 1
            else
                call type2_sph%move(type1_sph, mix, mix_ePot, type2_obs%Nrej)
                type2_obs%Nmove = type2_obs%Nmove + 1
            end if            
            
        end do
        
        call type1_sph%widom(type1_obs%activ)
        call type2_sph%widom(type2_obs%activ)
        
        call type1_obs%addReject()
        call type2_obs%addReject()
        
        if (iStep <= Ntherm) then
        
            call type1_sph%adaptDx(iStep, type1_obs%rejRateSum, &
                type1_io%report)
            write(type1_io%dx, *) iStep, type1_sph%getDx(), &
                type1_obs%rejRateSum/real(iStep, DP)
            write(type1_io%obsTherm, *) iStep, type1_obs%ePot, &
                type1_obs%activ
        
            call type2_sph%adaptDx(iStep, type2_obs%rejRateSum, &
                type2_io%report)
            write(type2_io%dx, *) iStep, type2_sph%getDx(), &
                type2_obs%rejRateSum/real(iStep, DP)
            write(type2_io%obsTherm, *) iStep, type2_obs%ePot, &
                type2_obs%activ
                
            write(unitObsTherm, *) iStep, type1_obs%ePot + type2_obs%ePot + &
                mix_ePot
            write(mix_unitObsTherm, *) iStep, mix_ePot
        
        else
        
            call type1_obs%addPhysical()
            write(type1_io%obs, *) iStep, type1_obs%ePot, &
                type1_obs%activ
        
            call type2_obs%addPhysical()
            write(type2_io%obs, *) iStep, type2_obs%ePot, &
                type2_obs%activ                
                
            mix_ePotSum = mix_ePotSum + mix_ePot            
            write(unitObs, *) iStep, type1_obs%ePot + type2_obs%ePot + mix_ePot
            write(mix_unitObs, *) iStep, mix_ePot

            if (snap) then
                call type1_sph%snapShot(type1_io%snapShots)
                call type2_sph%snapShot(type2_io%snapShots)
            end if
            
        end if
    
    end do
    call cpu_time(tFin)

    write(*, *) "End of cycles"

! End -----------------------------------------------------

    call type1_sph%overlapTest()
    call type1_sph%consistTest(type1_obs%ePot, type1_io%report)
    call type1_sph%snapShot(type1_io%snapFin)
    call type1_obs%results(type1_sph%getNcol(), type1_io%report)

    call type2_sph%overlapTest()
    call type2_sph%consistTest(type2_obs%ePot, type2_io%report)
    call type2_sph%snapShot(type2_io%snapFin)
    call type2_obs%results(type2_sph%getNcol(), type2_io%report)
    
    call mix%overlapTest(type1_sph%X, type2_sph%X)
    mix_ePot_total = mix%ePot_total(type1_sph%X, type2_sph%X)
    call mix_results(mix_ePot, mix_ePot_total, mix_ePotSum, mix_unitReport)
    
    ePot = type1_obs%ePot + type2_obs%ePot + mix_ePot
    ePot_total = type1_sph%ePot_total() + type2_sph%ePot_total() + &
        mix_ePot_total
    ePotSum = type1_obs%ePotSum + type2_obs%ePotSum + mix_ePotSum
    call results(ePot, ePot_total, ePotSum, tFin-tIni, unitReport)
    
    close(unitReport)
    close(unitObsTherm)
    close(unitObs)
    
    call type1_sph%destroy()
    call type1_io%close()
    
    call type2_sph%destroy()
    call type2_io%close()
    
    call mix%destroy()
    close(mix_unitReport)
    close(mix_unitObsTherm)
    close(mix_unitObs)
    
end program mc_canonical
