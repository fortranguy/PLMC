!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program mc_canonical

use, intrinsic :: iso_fortran_env
use data_constants
use data_mc
use data_distrib
use class_mixingPotential
use class_dipolarSpheres
use class_hardSpheres
use class_observables
use class_units
use mod_tools

implicit none
    
    ! Monte-Carlo variables
    integer :: iStep, iMove, iRotate !< Monte-Carlo counters
    integer :: iColRand !< random particle
    real(DP) :: rand !< random number in between 0 and 1
    real(DP) :: tIni, tFin !< CPU initial and final time
    
    ! Total physical system variables
    real(DP) :: Epot, EpotSum !< potential energy : at a Monte-Carlo step
    real(DP) :: Epot_conf !< potential energy : complete calculation from a configuration
    integer :: report_unit  !< data & results file
    integer :: obsThermal_unit !< observables files : Thermalisation
    integer :: obsEquilib_unit !< observables files : Equilibrium
    
    ! Mixing potential between 2 types
    type(MixingPotential) :: mix !< short-range potential
    real(DP) :: mix_Epot, mix_EpotSum
    real(DP) :: mix_Epot_conf
    integer :: mix_report_unit
    integer :: mix_Epot_unit !< tabulated potential file
    integer :: mix_obsThermal_unit
    integer :: mix_obsEquilib_unit
    
    ! Type 1 : Dipolar spheres : Ewald summation
    type(DipolarSpheres) :: type1_sph !< physical properties and Monte-Carlo subroutines
    type(MoreObservables) :: type1_obs !< energy & inverse of activity (-> chemical potential)
    type(MoreUnits) :: type1_io !< (input/)output files
    
    ! Type 2 : Hard spheres
    type(HardSpheres) :: type2_sph
    type(Observables) :: type2_obs
    type(Units) :: type2_io
    
! Beginning ----------------------------------------------------------------------------------------
    
    write(output_unit, *) "Monte-Carlo Mix - Canonical : Volume =", Volume

    ! Initialisations & reports
    
    open(newunit=report_unit, recl=4096, file="report.txt", status='new', action='write')
    open(newunit=obsThermal_unit, recl=4096, file="obsThermal.out", status='new', &
         action='write')
    open(newunit=obsEquilib_unit, recl=4096, file="obsEquilib.out", status='new', &
         action='write')
    call report(report_unit)
    !call initRandomSeed(report_unit)
    
    call mix%construct()
    mix_EpotSum = 0._DP
    open(newunit=mix_report_unit, recl=4096, file="mix_report.txt", status='new', action='write')
    open(newunit=mix_Epot_unit, recl=4096, file="mix_Epot.out", status='new', action='write')
    open(newunit=mix_obsThermal_unit, recl=4096, file="mix_obsThermal.out", &
         status='new', action='write')
    open(newunit=mix_obsEquilib_unit, recl=4096, file="mix_obsEquilib.out", status='new', &
         action='write')
    call mix%report(mix_report_unit)
    call mix%Epot_print(mix_Epot_unit)
    
    call type1_sph%construct(mix%getCell_Lsize(), mix%getRcut()) !< type1_sph needs mix%rCut for the 
                                                                 !< Cell List method
    call type1_obs%init()
    call type1_io%open(type1_sph%getName())
    call type1_sph%report(type1_io%report)
    call type1_sph%printInfo(type1_io%report)
    call type1_sph%Epot_real_print(type1_io%Epot)
    
    call type2_sph%construct(mix%getCell_Lsize(), mix%getRcut())
    call type2_obs%init()
    call type2_io%open(type2_sph%getName())
    call type2_sph%report(type2_io%report)
    call type2_sph%printInfo(type2_io%report)
    call type2_sph%Epot_print(type2_io%Epot)
    
    ! Initial condition
    
    call initialCondition(type1_sph, type2_sph, mix%getRmin(), report_unit)
    
    call type1_sph%overlapTest()
    call type1_sph%Epot_reci_init()
    type1_obs%Epot = type1_sph%Epot_conf()
    call type1_sph%snapShot_positions(0, type1_io%snapIni_positions)
    call type1_sph%snapShot_orientations(0, type1_io%snapIni_orientations)
    call type1_sph%cols_to_cells(type2_sph%positions) !< Cell List : filling cells with particles
    
    call type2_sph%overlapTest()
    type2_obs%Epot = type2_sph%Epot_conf()
    call type2_sph%snapShot_positions(0, type2_io%snapIni_positions)
    call type2_sph%cols_to_cells(type1_sph%positions)
    
    call mix%overlapTest(type1_sph%positions, type2_sph%positions)
    mix_Epot = mix%Epot_conf(type1_sph%positions, type2_sph%positions)
    
    Epot_conf = type1_obs%Epot + type2_obs%Epot + mix_Epot
    write(output_unit, *) "Initial potential energy =", Epot_conf
    write(obsThermal_unit, *) 0, Epot_conf
    
! Middle -------------------------------------------------------------------------------------------
        
    write(output_unit, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    MC_Cycle : do iStep = 1, Ntherm + Nstep
    
        ! Ewald summation : reinitialize the structure factor to prevent it from drifting.
        if (modulo(iStep, type1_sph%getStructure_iStep()) == 0) then
            call type1_sph%Epot_reci_structure_reInit(iStep, type1_io%structure_moduli)
        end if
    
        MC_Move : do iMove = 1, Nmove
        
            ! Randomly choosing a particle among both types
            call random_number(rand)
            iColRand = int(rand*real(Ncol, DP)) + 1
            
            ! Moving a particle : 
            if (iColRand <= type1_sph%getNcol()) then
                call type1_sph%move(iColRand, type2_sph, mix, type1_obs%Epot, mix_Epot, &
                                    type1_obs%Nrej)
                type1_obs%Nmove = type1_obs%Nmove + 1
            else
                iColRand = iColRand - type1_sph%getNcol()
                call type2_sph%move(iColRand, type1_sph, mix, type2_obs%Epot, mix_Epot, &
                                    type2_obs%Nrej)
                type2_obs%Nmove = type2_obs%Nmove + 1
            end if
            
        end do MC_Move
        
        MC_Rotate : do iRotate = 1, Nrotate
        
            call random_number(rand)
            iColRand = int(rand*real(type1_sph%getNcol(), DP)) + 1
 
            call type1_sph%rotate(iColRand, type1_obs%Epot, type1_obs%NrejRot)
            type1_obs%Nrotate = type1_obs%Nrotate + 1
            
        end do MC_Rotate
        
        ! Rejections rates updates
        type1_obs%rej = real(type1_obs%Nrej, DP)/real(type1_obs%Nmove, DP)
        type1_obs%Nrej = 0; type1_obs%Nmove = 0
        
        type1_obs%rejRot = real(type1_obs%NrejRot, DP)/real(type1_obs%Nrotate, DP)
        type1_obs%NrejRot = 0; type1_obs%Nrotate = 0

        type2_obs%rej = real(type2_obs%Nrej, DP)/real(type2_obs%Nmove, DP)
        type2_obs%Nrej = 0; type2_obs%Nmove = 0
        
        MC_Regime : if (iStep <= Ntherm) then ! Thermalisation
        
            ! Initial displacements & rejections
            if (iStep == 1) then
                write(type1_io%deltaX, *) iStep, type1_sph%getDeltaX(), type1_obs%rej
                write(type1_io%deltaM, *) iStep, type1_sph%getDeltaM(), type1_obs%rejRot
                write(type2_io%deltaX, *) iStep, type2_sph%getDeltaX(), type2_obs%rej
            end if
            
            ! Displacements adaptation           
            if (mod(iStep, type1_sph%getNadapt()) /= 0) then ! Rejections accumulation
                type1_obs%rejAdapt = type1_obs%rejAdapt + type1_obs%rej
            else ! Average & adaptation
                type1_obs%rejAdapt = type1_obs%rejAdapt/real(type1_sph%getNadapt()-1)
                call type1_sph%adaptDeltaX(type1_obs%rejAdapt)
                write(type1_io%deltaX, *) iStep, type1_sph%getDeltaX(), type1_obs%rejAdapt
                type1_obs%rejAdapt = 0._DP
            end if
            
            ! Rotation adaptation
            if (mod(iStep, type1_sph%getNadaptRot()) /= 0) then
                type1_obs%rejRotAdapt = type1_obs%rejRotAdapt + type1_obs%rejRot
            else
                type1_obs%rejRotAdapt = type1_obs%rejRotAdapt/real(type1_sph%getNadaptRot()-1)
                call type1_sph%adaptDeltaM(type1_obs%rejRotAdapt)
                write(type1_io%deltaM, *) iStep, type1_sph%getDeltaM(), type1_obs%rejRotAdapt
                type1_obs%rejRotAdapt = 0._DP
            end if

            if (mod(iStep, type2_sph%getNadapt()) /= 0) then
                type2_obs%rejAdapt = type2_obs%rejAdapt + type2_obs%rej
            else                
                type2_obs%rejAdapt = type2_obs%rejAdapt/real(type2_sph%getNadapt()-1)
                call type2_sph%adaptDeltaX(type2_obs%rejAdapt)
                write(type2_io%deltaX, *) iStep, type2_sph%getDeltaX(), type2_obs%rejAdapt
                type2_obs%rejAdapt = 0._DP
            end if
            
            ! Observables writing
            write(type1_io%obsThermal, *) iStep, type1_obs%Epot, 0._DP, type1_obs%rej, &
                                                        type1_obs%rejRot
            write(type2_io%obsThermal, *) iStep, type2_obs%Epot, 0._DP, type2_obs%rej
            write(mix_obsThermal_unit, *) iStep, mix_Epot
            write(obsThermal_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot
            
            if (iStep == Ntherm) then ! Definite thermalised displacements
                call type1_sph%definiteDeltaX(type1_obs%rej, type1_io%report)
                call type1_sph%definiteDeltaM(type1_obs%rejRot, type1_io%report)
                call type2_sph%definiteDeltaX(type2_obs%rej, type2_io%report)
            end if       
        
        else MC_Regime ! Thermalisation over -> Equilibrium
        
            ! Chemical potentials : Widom method
            call type1_sph%widom(type2_sph%positions, mix, type1_obs%activ)
            call type2_sph%widom(type1_sph%positions, mix, type2_obs%activ)
        
            ! Observables accumulations
            type1_obs%EpotSum = type1_obs%EpotSum + type1_obs%Epot
            type1_obs%activSum = type1_obs%activSum + type1_obs%activ
            type1_obs%rejSum = type1_obs%rejSum + type1_obs%rej
            type1_obs%rejRotSum = type1_obs%rejRotSum + type1_obs%rejRot
        
            type2_obs%EpotSum = type2_obs%EpotSum + type2_obs%Epot
            type2_obs%activSum = type2_obs%activSum + type2_obs%activ
            type2_obs%rejSum = type2_obs%rejSum + type2_obs%rej
                
            mix_EpotSum = mix_EpotSum + mix_Epot

            write(type1_io%obsEquilib, *) iStep, type1_obs%Epot, type1_obs%activ, type1_obs%rej, &
                                                     type1_obs%rejRot
            write(type2_io%obsEquilib, *) iStep, type2_obs%Epot, type2_obs%activ, type2_obs%rej
            write(mix_obsEquilib_unit, *) iStep, mix_Epot
            write(obsEquilib_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot

            if (snap) then ! Snap shots of the configuration
                call type1_sph%snapShot_positions(iStep, type1_io%snapShots_positions)
                call type1_sph%snapShot_orientations(iStep, type1_io%snapShots_orientations)
                call type2_sph%snapShot_positions(iStep, type2_io%snapShots_positions)
            end if
            
        end if MC_Regime
    
    end do MC_Cycle
    call cpu_time(tFin)

    write(output_unit, *) "End of cycles"

! End ----------------------------------------------------------------------------------------------

    ! Tests & results

    call type1_sph%overlapTest()
    call type1_sph%Epot_reci_init()
    call type1_sph%consistTest(type1_obs%Epot, type1_io%report)
    call type1_sph%snapShot_positions(0, type1_io%snapFin_positions)
    call type1_sph%snapShot_orientations(0, type1_io%snapFin_orientations)
    call type1_obs%results(type1_sph%getNcol(), type1_io%report)
    
    call type2_sph%overlapTest()
    call type2_sph%consistTest(type2_obs%Epot, type2_io%report)
    call type2_sph%snapShot_positions(0, type2_io%snapFin_positions)
    call type2_obs%results(type2_sph%getNcol(), type2_io%report)
    
    call mix%overlapTest(type1_sph%positions, type2_sph%positions)
    mix_Epot_conf = mix%Epot_conf(type1_sph%positions, type2_sph%positions)
    call consistTest(mix_Epot, mix_Epot_conf, mix_report_unit)
    call mix_results(mix_EpotSum, mix_report_unit)
    
    Epot = type1_obs%Epot + type2_obs%Epot + mix_Epot
    Epot_conf = type1_sph%Epot_conf() + type2_sph%Epot_conf() + mix_Epot_conf
    call consistTest(Epot, Epot_conf, report_unit)
    EpotSum = type1_obs%EpotSum + type2_obs%EpotSum + mix_EpotSum
    call results(EpotSum, tFin-tIni, report_unit)
    
    ! Finalisations
    
    call type1_sph%destroy()
    call type1_io%close()
    
    call type2_sph%destroy()
    call type2_io%close()
    
    call mix%destroy()
    close(mix_report_unit)
    close(mix_Epot_unit)
    close(mix_obsThermal_unit)
    close(mix_obsEquilib_unit)
    
    close(report_unit)
    close(obsThermal_unit)
    close(obsEquilib_unit)
    
end program mc_canonical
