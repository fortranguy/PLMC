!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use iso_fortran_env
use data_constants
use data_mc
use data_distrib
use class_mixingPotential
use class_interactingSpheres
use class_hardSpheres
use class_observables
use class_units
use mod_tools

implicit none

! Beginning ----------------------------------------------------------------------------------------

    ! Declarations
    
    integer :: iStep, iMove !< Monte-Carlo counters
    integer :: iColRand !< random choice of a particle
    real(DP) :: rand !< random number in between 0 and 1
    real(DP) :: tIni, tFin !< initial and final time
    
    !   Type 1 : Interacting spheres
    type(InteractingSpheres) :: type1_sph !< Monte-Carlo subroutines
    type(Observables) :: type1_obs !< e.g. Energy
    type(Units) :: type1_io        !< input/output files
    
    !   Type 2 : Hard spheres
    type(HardSpheres) :: type2_sph
    type(Observables) :: type2_obs
    type(Units) :: type2_io
    
    !   Mixing between 2 types
    type(MixingPotential) :: mix
    real(DP) :: mix_Epot, mix_EpotSum, mix_Epot_conf !< potential energy
    integer, parameter :: mix_report_unit = 10 !< data & results
    integer, parameter :: mix_obsTherm_unit = 11, mix_obs_unit = 12 !< observable(s)
    
    !   Whole system variables
    real(DP) :: Epot, EpotSum, Epot_conf
    integer, parameter :: report_unit = 13
    integer, parameter :: obsTherm_unit = 14, obs_unit = 15
    
    write(output_unit, *) "Monte-Carlo Mix - Canonical : Volume =", product(Lsize)

    ! Initialisations & reports
    
    call type1_sph%construct()
    call type1_obs%init()
    call type1_io%open(type1_sph%getName())
    
    call type2_sph%construct()
    call type2_obs%init()
    call type2_io%open(type2_sph%getName())
    
    call mix%construct()
    mix_EpotSum = 0._DP
    open(unit=mix_report_unit, recl=4096, file="mix_report.out", status='new', action='write')
    open(unit=mix_obsTherm_unit, recl=4096, file="mix_obsTherm.out", status='new', action='write')
    open(unit=mix_obs_unit, recl=4096, file="mix_obs.out", status='new', action='write')
        
    open(unit=report_unit, recl=4096, file="report.out", status='new', action='write')
    open(unit=obsTherm_unit, recl=4096, file="obsTherm.out", status='new', action='write')
    open(unit=obs_unit, recl=4096, file="obs.out", status='new', action='write')
    
    call type1_sph%report(type1_io%report)
    call type1_sph%printInfo(type1_io%report)    
    call type2_sph%report(type2_io%report)
    call type2_sph%printInfo(type2_io%report)
    
    call mix%report(mix_report_unit)
    call report(report_unit)
    
    call initRandomSeed(report_unit)
    
    ! Initial condition
    
    call initialCondition(type1_sph, type2_sph, mix%getRmin(), report_unit)
    
    call type1_sph%overlapTest()
    call type2_sph%overlapTest()
    call mix%overlapTest(type1_sph%X, type2_sph%X)
    
    type1_obs%Epot = type1_sph%Epot_conf()
    call type1_sph%snapShot(type1_io%snapIni)
    call type1_sph%cols_to_cells(type2_sph%X)    
    
    type2_obs%Epot = type2_sph%Epot_conf()
    call type2_sph%snapShot(type2_io%snapIni)
    call type2_sph%cols_to_cells(type1_sph%X)    
    
    mix_Epot = mix%Epot_conf(type1_sph%X, type2_sph%X)    
    Epot_conf = type1_obs%Epot + type2_obs%Epot + mix_Epot
    
! Middle -------------------------------------------------------------------------------------------
        
    write(output_unit, *) "Beginning of cycles"    
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        ! Moving a particle : Metropolis algorithm
        do iMove = 1, Nmove
        
            call random_number(rand)
            iColRand = int(rand*real(Ncol, DP)) + 1
            if (iColRand <= type1_sph%getNcol()) then
                call type1_sph%move(type2_sph, mix, type1_obs%Epot, mix_Epot, type1_obs%Nrej)
                type1_obs%Nmove = type1_obs%Nmove + 1
            else
                call type2_sph%move(type1_sph, mix, type2_obs%Epot, mix_Epot, type2_obs%Nrej)
                type2_obs%Nmove = type2_obs%Nmove + 1
            end if
            
        end do
        
        ! Widom method : chemical potentials
        call type1_sph%widom(type1_obs%activ)
        call type2_sph%widom(type2_obs%activ)
        
        ! Rejections rate
        type1_obs%rej = real(type1_obs%Nrej, DP)/real(type1_obs%Nmove, DP)
        type1_obs%Nrej = 0; type1_obs%Nmove = 0        
        type2_obs%rej = real(type2_obs%Nrej, DP)/real(type2_obs%Nmove, DP)
        type2_obs%Nrej = 0; type2_obs%Nmove = 0
        
        if (iStep <= Ntherm) then ! Thermalisation
        
            if (iStep==1) then
                write(type1_io%dx, *) iStep, type1_sph%getDx(), type1_obs%rej
                write(type2_io%dx, *) iStep, type2_sph%getDx(), type2_obs%rej
            end if
            
            ! Displacements optimisation : type 1            
            if (mod(iStep, type1_sph%getNadapt()) == 0) then
                type1_obs%rejAdapt = type1_obs%rejAdapt/real(type1_sph%getNadapt()-1)
                call type1_sph%adaptDx(type1_obs%rejAdapt)
                write(type1_io%dx, *) iStep, type1_sph%getDx(), type1_obs%rejAdapt
                type1_obs%rejAdapt = 0._DP
            else
                type1_obs%rejAdapt = type1_obs%rejAdapt + type1_obs%rej
            end if
            ! Displacements optimisation : type 2
            if (mod(iStep, type2_sph%getNadapt()) == 0) then
                type2_obs%rejAdapt = type2_obs%rejAdapt/real(type2_sph%getNadapt()-1)
                call type2_sph%adaptDx(type2_obs%rejAdapt)
                write(type2_io%dx, *) iStep, type2_sph%getDx(), type2_obs%rejAdapt
                type2_obs%rejAdapt = 0._DP
            else
                type2_obs%rejAdapt = type2_obs%rejAdapt + type2_obs%rej
            end if
            
            ! Observables writing
            write(type1_io%obsTherm, *) iStep, type1_obs%Epot, type1_obs%activ, type1_obs%rej
            write(type2_io%obsTherm, *) iStep, type2_obs%Epot, type2_obs%activ, type2_obs%rej
            write(mix_obsTherm_unit, *) iStep, mix_Epot
            write(obsTherm_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot
            
            if (iStep == Ntherm) then ! Definite thermalised displacement
                call type1_sph%definiteDx(type1_obs%rej, type1_io%report)
                call type2_sph%definiteDx(type2_obs%rej, type2_io%report)
            end if       
        
        else ! Observables accumulations & writing
        
            type1_obs%EpotSum = type1_obs%EpotSum + type1_obs%Epot
            type1_obs%activSum = type1_obs%activSum + type1_obs%activ
            type1_obs%rejSum = type1_obs%rejSum + type1_obs%rej
            write(type1_io%obs, *) iStep, type1_obs%Epot, type1_obs%activ, type1_obs%rej
        
            type2_obs%EpotSum = type2_obs%EpotSum + type2_obs%Epot
            type2_obs%activSum = type2_obs%activSum + type2_obs%activ
            type2_obs%rejSum = type2_obs%rejSum + type2_obs%rej
            write(type2_io%obs, *) iStep, type2_obs%Epot, type2_obs%activ, type2_obs%rej
                
            mix_EpotSum = mix_EpotSum + mix_Epot
            write(mix_obs_unit, *) iStep, mix_Epot
            write(obs_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot

            if (snap) then ! snap shots of the configuration
                call type1_sph%snapShot(type1_io%snapShots)
                call type2_sph%snapShot(type2_io%snapShots)
            end if
            
        end if
    
    end do
    call cpu_time(tFin)

    write(output_unit, *) "End of cycles"

! End ----------------------------------------------------------------------------------------------

    ! Tests & results

    call type1_sph%overlapTest()
    call type1_sph%consistTest(type1_obs%Epot, type1_io%report)
    call type1_sph%snapShot(type1_io%snapFin)
    call type1_obs%results(type1_sph%getNcol(), type1_io%report)

    call type2_sph%overlapTest()
    call type2_sph%consistTest(type2_obs%Epot, type2_io%report)
    call type2_sph%snapShot(type2_io%snapFin)
    call type2_obs%results(type2_sph%getNcol(), type2_io%report)
    
    call mix%overlapTest(type1_sph%X, type2_sph%X)
    mix_Epot_conf = mix%Epot_conf(type1_sph%X, type2_sph%X)
    call mix_results(mix_Epot, mix_Epot_conf, mix_EpotSum, mix_report_unit)
    
    Epot = type1_obs%Epot + type2_obs%Epot + mix_Epot
    Epot_conf = type1_sph%Epot_conf() + type2_sph%Epot_conf() + mix_Epot_conf
    EpotSum = type1_obs%EpotSum + type2_obs%EpotSum + mix_EpotSum
    call results(Epot, Epot_conf, EpotSum, tFin-tIni, report_unit)
    
    ! Finalisations
    
    call type1_sph%destroy()
    call type1_io%close()    
    call type2_sph%destroy()
    call type2_io%close()
    
    call mix%destroy()
    close(mix_report_unit)
    close(mix_obsTherm_unit)
    close(mix_obs_unit)
    
    close(report_unit)
    close(obsTherm_unit)
    close(obs_unit)
    
end program mc_canonical
