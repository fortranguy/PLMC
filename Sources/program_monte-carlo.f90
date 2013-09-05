!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program monteCarlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_particles, only : Ncol
use data_monteCarlo, only : Nstep, Ntherm, Nmove, Nrotate
use data_distribution, only : snap
use class_observables
use class_mixingPotential
use class_dipolarSpheres
use class_hardSpheres
use class_units
use module_tools, only : initRandomSeed, initialCondition, report, consistTest, printResults, &
                         mix_printResults

implicit none
    
    ! Monte-Carlo variables
    integer :: iStep !< step counter
    integer :: iChangeRand !< random change
    integer :: iChange !< change counters
    integer :: iColRand !< random particle
    real(DP) :: random !< random number between 0 and 1
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
    type(DipolarSpheres) :: type1_spheres !< physical properties and Monte-Carlo subroutines
    type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
    type(MoreUnits) :: type1_units !< files units
    
    ! Type 2 : Hard spheres
    type(HardSpheres) :: type2_spheres
    type(Observables) :: type2_obs
    type(Units) :: type2_units
    
! Beginning ----------------------------------------------------------------------------------------
    
    write(output_unit, *) "Monte-Carlo Simulation : Canonical ensemble"

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
    open(newunit=mix_Epot_unit, recl=4096, file="mix_Epot.tmp", status='new', action='write')
    open(newunit=mix_obsThermal_unit, recl=4096, file="mix_obsThermal.out", &
         status='new', action='write')
    open(newunit=mix_obsEquilib_unit, recl=4096, file="mix_obsEquilib.out", status='new', &
         action='write')
    call mix%Epot_print(mix_Epot_unit)
    call mix%printReport(mix_report_unit)
    
    call type1_spheres%construct(mix%getCell_size(), mix%getRcut())
    call type1_obs%init()
    call type1_units%open(type1_spheres%getName())
    call type1_spheres%Epot_real_print(type1_units%Epot)
    call type1_spheres%Epot_reci_countNwaveVectors(type1_units%waveVectors)
    call type1_spheres%printDensity(type1_units%report)
    call type1_spheres%printReport(type1_units%report)
    
    call type2_spheres%construct(mix%getCell_size(), mix%getRcut())
    call type2_obs%init()
    call type2_units%open(type2_spheres%getName())
    call type2_spheres%Epot_print(type2_units%Epot)
    call type2_spheres%printDensity(type2_units%report)
    call type2_spheres%printReport(type2_units%report)
    
    ! Initial condition
    
    call initialCondition(type1_spheres, type2_spheres, mix%getRmin(), report_unit)
    
    call type1_spheres%overlapTest()
    call type1_spheres%Epot_reci_init()
    call type1_spheres%Epot_bound_totalMoment_init()
    type1_obs%Epot = type1_spheres%Epot_conf()
    call type1_spheres%snapShot_positions(0, type1_units%snapIni_positions)
    call type1_spheres%snapShot_orientations(0, type1_units%snapIni_orientations)
    call type1_spheres%all_cols_to_cells(type2_spheres) !< filling cells with particles
    
    call type2_spheres%overlapTest()
    type2_obs%Epot = type2_spheres%Epot_conf()
    call type2_spheres%snapShot_positions(0, type2_units%snapIni_positions)
    call type2_spheres%all_cols_to_cells(type1_spheres)
    
    call mix%overlapTest(type1_spheres, type2_spheres)
    mix_Epot = mix%Epot_conf(type1_spheres, type2_spheres)
    
    Epot_conf = type1_obs%Epot + type2_obs%Epot + mix_Epot
    write(output_unit, *) "Initial potential energy =", Epot_conf
    write(obsThermal_unit, *) 0, Epot_conf
    
! Middle -------------------------------------------------------------------------------------------
        
    write(output_unit, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    MC_Cycle : do iStep = 1, Ntherm + Nstep
        
        MC_Change : do iChange = 1, Nmove + Nrotate
        
            ! Randomly choosing the change
            call random_number(random)
            iChangeRand = int(random*real(Nmove+Nrotate, DP)) + 1
            
            if (iChangeRand <= Nmove) then ! change = move

                ! Randomly choosing a particle among both types
                call random_number(random)
                iColRand = int(random*real(Ncol, DP)) + 1
                
                ! Moving a particle : 
                if (iColRand <= type1_spheres%getNcol()) then
                    call type1_spheres%move(iColRand, type2_spheres, mix, type1_obs, mix_Epot)
                    type1_obs%Nmove = type1_obs%Nmove + 1
                else
                    iColRand = iColRand - type1_spheres%getNcol()
                    call type2_spheres%move(iColRand, type1_spheres, mix, type2_obs, mix_Epot)
                    type2_obs%Nmove = type2_obs%Nmove + 1
                end if
                
            else ! change = rotate
                
                call random_number(random)
                iColRand = int(random*real(type1_spheres%getNcol(), DP)) + 1
     
                call type1_spheres%rotate(iColRand, type1_obs)
                type1_obs%Nrotate = type1_obs%Nrotate + 1
                
            end if
            
        end do MC_Change
        
        ! Rejections rates updates
        type1_obs%move_reject = real(type1_obs%move_Nreject, DP)/real(type1_obs%Nmove, DP)
        type1_obs%move_Nreject = 0; type1_obs%Nmove = 0
        
        type1_obs%rotate_reject = real(type1_obs%rotate_Nreject, DP)/real(type1_obs%Nrotate, DP)
        type1_obs%rotate_Nreject = 0; type1_obs%Nrotate = 0

        type2_obs%move_reject = real(type2_obs%move_Nreject, DP)/real(type2_obs%Nmove, DP)
        type2_obs%move_Nreject = 0; type2_obs%Nmove = 0
        
        MC_Regime : if (iStep <= Ntherm) then ! Thermalisation
            
            ! Displacements adaptation
            if (mod(iStep, type1_spheres%getMove_Nadapt()) /= 0) then ! Rejections accumulation
                type1_obs%move_rejectAdapt = type1_obs%move_rejectAdapt + type1_obs%move_reject
            else ! Average & adaptation
                type1_obs%move_rejectAdapt = type1_obs%move_rejectAdapt / &
                                             real(type1_spheres%getMove_Nadapt()-1)
                call type1_spheres%adaptDeltaX(type1_obs%move_rejectAdapt)
                write(type1_units%move_delta, *) iStep, type1_spheres%getMove_delta(), &
                                             type1_obs%move_rejectAdapt
                type1_obs%move_rejectAdapt = 0._DP
            end if
            
            ! Rotation adaptation
            if (mod(iStep, type1_spheres%getRotate_Nadapt()) /= 0) then
                type1_obs%rotate_rejectAdapt = type1_obs%rotate_rejectAdapt + type1_obs%rotate_reject
            else
                type1_obs%rotate_rejectAdapt = type1_obs%rotate_rejectAdapt / &
                                               real(type1_spheres%getRotate_Nadapt()-1)
                call type1_spheres%adaptDeltaM(type1_obs%rotate_rejectAdapt)
                write(type1_units%deltaM, *) iStep, type1_spheres%getDeltaM(), &
                                             type1_obs%rotate_rejectAdapt
                type1_obs%rotate_rejectAdapt = 0._DP
            end if

            if (mod(iStep, type2_spheres%getMove_Nadapt()) /= 0) then
                type2_obs%move_rejectAdapt = type2_obs%move_rejectAdapt + type2_obs%move_reject
            else                
                type2_obs%move_rejectAdapt = type2_obs%move_rejectAdapt / &
                                             real(type2_spheres%getMove_Nadapt()-1)
                call type2_spheres%adaptDeltaX(type2_obs%move_rejectAdapt)
                write(type2_units%move_delta, *) iStep, type2_spheres%getMove_delta(), &
                                             type2_obs%move_rejectAdapt
                type2_obs%move_rejectAdapt = 0._DP
            end if
            
            ! Observables writing
            write(type1_units%obsThermal, *) iStep, type1_obs%Epot, 0._DP, type1_obs%move_reject, &
                                             type1_obs%rotate_reject
            write(type2_units%obsThermal, *) iStep, type2_obs%Epot, 0._DP, type2_obs%move_reject
            write(mix_obsThermal_unit, *) iStep, mix_Epot
            write(obsThermal_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot
            
            if (iStep == Ntherm) then ! Definite thermalised displacements
                call type1_spheres%definiteDeltaX(type1_obs%move_reject, type1_units%report)
                call type1_spheres%definiteDeltaM(type1_obs%rotate_reject, type1_units%report)
                call type2_spheres%definiteDeltaX(type2_obs%move_reject, type2_units%report)
            end if       
        
        else MC_Regime ! Thermalisation over -> Equilibrium
        
            ! Chemical potentials : Widom method
            call type1_spheres%widom(type2_spheres%positions, mix, type1_obs%activ)
            call type2_spheres%widom(type1_spheres%positions, mix, type2_obs%activ)
        
            ! Observables accumulations
            type1_obs%EpotSum = type1_obs%EpotSum + type1_obs%Epot
            type1_obs%activSum = type1_obs%activSum + type1_obs%activ
            type1_obs%move_rejectSum = type1_obs%move_rejectSum + type1_obs%move_reject
            type1_obs%rotate_rejectSum = type1_obs%rotate_rejectSum + type1_obs%rotate_reject
        
            type2_obs%EpotSum = type2_obs%EpotSum + type2_obs%Epot
            type2_obs%activSum = type2_obs%activSum + type2_obs%activ
            type2_obs%move_rejectSum = type2_obs%move_rejectSum + type2_obs%move_reject
                
            mix_EpotSum = mix_EpotSum + mix_Epot

            write(type1_units%obsEquilib, *) iStep, type1_obs%Epot, type1_obs%activ, &
                                             type1_obs%move_reject, type1_obs%rotate_reject
            write(type2_units%obsEquilib, *) iStep, type2_obs%Epot, type2_obs%activ, &
                                             type2_obs%move_reject
            write(mix_obsEquilib_unit, *) iStep, mix_Epot
            write(obsEquilib_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot

            if (snap) then ! Snap shots of the configuration
                call type1_spheres%snapShot_positions(iStep, type1_units%snapShots_positions)
                call type1_spheres%snapShot_orientations(iStep, type1_units%snapShots_orientations)
                call type2_spheres%snapShot_positions(iStep, type2_units%snapShots_positions)
            end if
            
        end if MC_Regime
        
        ! Ewald summation : reinitialize the structure factor to prevent it from drifting.
        if (modulo(iStep, type1_spheres%getStructure_iStep()) == 0) then
            call type1_spheres%Epot_reci_structure_reInit(iStep, type1_units%structure_modulus)
        end if
        
        ! Boundary conditions : reinitialize the total moment to prevent it from drifting.
        if (modulo(iStep, type1_spheres%getTotalMoment_iStep()) == 0) then
            call type1_spheres%Epot_bound_totalMoment_reInit(iStep, type1_units%totalMoment_modulus)
        end if
    
    end do MC_Cycle
    call cpu_time(tFin)

    write(output_unit, *) "End of cycles"

! End ----------------------------------------------------------------------------------------------

    ! Tests & results

    call type1_spheres%overlapTest()
    call type1_spheres%Epot_reci_init()
    call type1_spheres%Epot_bound_totalMoment_init()
    call type1_spheres%consistTest(type1_obs%Epot, type1_units%report)
    call type1_spheres%snapShot_positions(0, type1_units%snapFin_positions)
    call type1_spheres%snapShot_orientations(0, type1_units%snapFin_orientations)
    call type1_obs%printResults(type1_spheres%getNcol(), type1_units%report)
    
    call type2_spheres%overlapTest()
    call type2_spheres%consistTest(type2_obs%Epot, type2_units%report)
    call type2_spheres%snapShot_positions(0, type2_units%snapFin_positions)
    call type2_obs%printResults(type2_spheres%getNcol(), type2_units%report)
    
    call mix%overlapTest(type1_spheres, type2_spheres)
    mix_Epot_conf = mix%Epot_conf(type1_spheres, type2_spheres)
    call consistTest(mix_Epot, mix_Epot_conf, mix_report_unit)
    call mix_printResults(mix_EpotSum, mix_report_unit)
    
    Epot = type1_obs%Epot + type2_obs%Epot + mix_Epot
    Epot_conf = type1_spheres%Epot_conf() + type2_spheres%Epot_conf() + mix_Epot_conf
    call consistTest(Epot, Epot_conf, report_unit)
    EpotSum = type1_obs%EpotSum + type2_obs%EpotSum + mix_EpotSum
    call printResults(EpotSum, tFin-tIni, report_unit)
    
    ! Finalisations
    
    call type1_spheres%destroy()
    call type1_units%close()
    
    call type2_spheres%destroy()
    call type2_units%close()
    
    call mix%destroy()
    close(mix_report_unit)
    close(mix_Epot_unit)
    close(mix_obsThermal_unit)
    close(mix_obsEquilib_unit)
    
    close(report_unit)
    close(obsThermal_unit)
    close(obsEquilib_unit)
    
end program monteCarlo_canonical_bulk
