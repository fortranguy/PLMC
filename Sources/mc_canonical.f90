!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use data_constants
use data_mc
use data_distrib
use mod_tools
use class_interactingSpheres
use class_hardSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none

! Beginning
    
    integer :: iStep, iMove
    integer :: iColRand
    real(DP) :: rand
    real(DP) :: tIni, tFin
    
    integer, parameter :: unitReport = 10
    
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

    ! Initialisation
    
    call type1_sph%construct()
    call type1_obs%init()
    call type1_io%open(type1_sph%getName())
    
    call type2_sph%construct()
    call type2_obs%init()
    call type2_io%open(type2_sph%getName())
    
    call mix%construct()

    write(*, *) "Monte-Carlo - Canonical : Volume =", product(Lsize)    
    
    open(unit=unitReport, recl=4096, file="report.out", status='new', &
    	action='write')
	call report(unitReport)
    call type1_sph%report(type1_io%report)
    call type2_sph%report(type2_io%report)
    
    call init_random_seed(unitReport)
    
    ! Initial condition
    
    call initialCondition(type1_sph%X, type1_sph%getRmin(), type1_io%report)
    call type1_sph%overlapTest()
    type1_obs%ePot_total = type1_sph%ePot_total()
    call type1_sph%snapShot(type1_io%snapIni)
    call type1_sph%cols_to_cells()
    
    call initialCondition(type2_sph%X, type2_sph%getRmin(), type2_io%report)
    call type2_sph%overlapTest()
    type2_obs%ePot_total = 0._DP
    call type2_sph%snapShot(type2_io%snapIni)
    call type2_sph%cols_to_cells()
    
! Middle --------------------------------------------------

    write(*, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    do iStep = 1, Ntherm + Nstep
    
        do iMove = 1, Nmove
        
            call random_number(rand)
            iColRand = int(rand*real(Ncol, DP)) + 1            
            if (iColRand <= type1_sph%getNcol()) then
                call type1_sph%move(type1_obs%ePot_total, type1_obs%Nrej)
                type1_obs%Nmove = type1_obs%Nmove + 1
            else
                call type2_sph%move(type2_obs%Nrej)
                type2_obs%Nmove = type2_obs%Nmove + 1
            end if            
            
        end do
        
        call type1_sph%widom(type1_obs%activExInv)
        call type2_sph%widom(type2_obs%activExInv)
        
        call type1_obs%addReject()
        call type2_obs%addReject()
        
        if (iStep <= Ntherm) then
        
            call type1_sph%adaptDx(iStep, type1_obs%rejRateSum, &
                type1_io%report)
            write(type1_io%dx, *) iStep, type1_sph%getDx(), &
                type1_obs%rejRateSum/real(iStep, DP)
            write(type1_io%obsTherm, *) iStep, type1_obs%ePot_total, &
                type1_obs%activExInv
        
            call type2_sph%adaptDx(iStep, type2_obs%rejRateSum, &
                type2_io%report)
            write(type2_io%dx, *) iStep, type2_sph%getDx(), &
                type2_obs%rejRateSum/real(iStep, DP)
            write(type2_io%obsTherm, *) iStep, type2_obs%ePot_total, &
                type2_obs%activExInv
        
        else
        
            call type1_obs%addPhysical()
            write(type1_io%obs, *) iStep, type1_obs%ePot_total, &
                type1_obs%activExInv
        
            call type2_obs%addPhysical()
            write(type2_io%obs, *) iStep, type2_obs%ePot_total, &
                type2_obs%activExInv

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
    call type1_sph%consistTest(type1_obs%ePot_total, type1_io%report)
    call type1_sph%snapShot(type1_io%snapFin)
    call type1_obs%results(tFin-tIni, type1_io%report)

    call type2_sph%overlapTest()
    call type2_sph%snapShot(type2_io%snapFin)
    call type2_obs%results(tFin-tIni, type2_io%report)

    close(unitReport)
    
    call type1_sph%destroy()
    call type1_io%close()
    
    call type2_sph%destroy()
    call type2_io%close()
    
    call mix%destroy()
    
end program mc_canonical
