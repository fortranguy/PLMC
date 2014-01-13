!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program monteCarlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_monteCarlo, only : decorrelFactor, Nthermal, Nadapt, Nstep
use data_distribution, only : snap
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units
use module_algorithms, only : move, widom, rotate
use module_tools, only : open_units, mix_open_units, init_randomSeed, set_initialCondition, &
                         print_report, mix_init, mix_final, init, final, adapt_move, adapt_rotate, &
                         test_consist, print_results

implicit none
    
    ! Monte-Carlo variables
    integer :: iStep, iChange, iChangeRand !< counters
    integer :: Ncol, iColRand !< number of particles & random particle
    integer :: Nmove, Nrotate !< number of changes
    real(DP) :: random !< random number between 0 and 1
    real(DP) :: tIni, tFin !< CPU initial and final time
    
    ! Total physical system variables
    real(DP) :: Epot, EpotSum !< potential energy : at a Monte-Carlo step
    real(DP) :: Epot_conf !< potential energy : complete calculation from a configuration
    integer :: report_unit !< data & results file
    integer :: obsThermal_unit, obsEquilib_unit !< observables files : thermalisation & equilibrium
    
    ! Mixing potential between 2 types
    type(MixingPotential) :: mix !< short-range potential
    real(DP) :: mix_Epot, mix_EpotSum, mix_Epot_conf
    integer :: mix_report_unit
    integer :: mix_Epot_unit !< tabulated potential file
    integer :: mix_obsThermal_unit, mix_obsEquilib_unit
    
    ! Type 1 : Dipolar spheres : Ewald summation
    type(DipolarSpheres) :: type1_spheres !< physical properties and Monte-Carlo subroutines
    type(MoreObservables) :: type1_obs !< e.g. energy, inverse of activity (-> chemical potential)
    type(MoreUnits) :: type1_units !< files units
    
    ! Type 2 : Hard spheres
    type(HardSpheres) :: type2_spheres
    type(Observables) :: type2_obs
    type(Units) :: type2_units

    call type1_spheres%construct()
    call type2_spheres%construct()
    call mix%construct(type1_spheres%get_sigma(), type2_spheres%get_sigma())
    
    Ncol = type1_spheres%get_Ncol() + type2_spheres%get_Ncol()
    Nmove = decorrelFactor * Ncol
    Nrotate = decorrelFactor * type1_spheres%get_Ncol()
    
! Beginning ----------------------------------------------------------------------------------------
    
    write(output_unit, *) "Monte-Carlo Simulation : Canonical ensemble"

    ! Initialisations & reports
    
    call open_units(report_unit, obsThermal_unit, obsEquilib_unit)
    call print_report(Ncol, Nmove, Nrotate, report_unit)
    !call init_randomSeed(report_unit)
    
    mix_EpotSum = 0._DP
    call mix_open_units(mix_report_unit, mix_Epot_unit, mix_obsThermal_unit, mix_obsEquilib_unit)
    
    call type1_obs%init()
    call type1_units%open(type1_spheres%get_name())
    call type1_spheres%print_density(Ncol, type1_units%report)
    
    call type2_obs%init()
    call type2_units%open(type2_spheres%get_name())
    call type2_spheres%print_density(Ncol, type2_units%report)
    
    ! Initial condition
    
    call set_initialCondition(type1_spheres, type2_spheres, mix%get_sigma(), report_unit)
    
    call mix_init(mix, type1_spheres, type2_spheres, mix_Epot_unit, mix_Epot)
    call mix%print_report(mix_report_unit)
    call init(type1_spheres, type2_spheres, mix, type1_units, type1_obs%Epot)
    call init(type2_spheres, type1_spheres, mix, type2_units, type2_obs%Epot)
    
    Epot_conf = type1_obs%Epot + type2_obs%Epot + mix_Epot
    write(output_unit, *) "Initial potential energy =", Epot_conf
    write(obsThermal_unit, *) 0, Epot_conf
    
! Middle -------------------------------------------------------------------------------------------
        
    write(output_unit, *) "Beginning of cycles"
    
    call cpu_time(tIni)
    MC_Cycle : do iStep = 1, Nthermal + Nstep
        
        MC_Change : do iChange = 1, Nmove + Nrotate
        
            ! Randomly choosing the change
            call random_number(random)
            iChangeRand = int(random*real(Nmove+Nrotate, DP)) + 1
            
            if (iChangeRand <= Nmove) then
            
                ! Randomly choosing the type
                call random_number(random)
                iColRand = int(random*real(Ncol, DP)) + 1
                
                if (iColRand <= type1_spheres%get_Ncol()) then
                    call move(type1_spheres, type1_obs, type2_spheres, mix, mix_Epot)
                else
                    call move(type2_spheres, type2_obs, type1_spheres, mix, mix_Epot)
                end if
                
            else
     
                call rotate(type1_spheres, type1_obs)
                
            end if
            
        end do MC_Change
        
        ! Rejections rates updates
        call type1_obs%update_rejections()
        call type2_obs%update_rejections()
        
        MC_Regime : if (iStep <= Nthermal) then ! Thermalisation
            
            ! Change adaptation
            if (mod(iStep, Nadapt) /= 0) then ! Rejections accumulation            
                type1_obs%move_rejectAdapt = type1_obs%move_rejectAdapt + type1_obs%move_reject
                type1_obs%rotate_rejectAdapt = type1_obs%rotate_rejectAdapt + type1_obs%rotate_reject
                type2_obs%move_rejectAdapt = type2_obs%move_rejectAdapt + type2_obs%move_reject                
            else ! Average & adaptation            
                call adapt_move(type1_spheres, iStep, type1_obs, type1_units%move_delta)
                call adapt_rotate(type1_spheres, iStep, type1_obs, type1_units%rotate_delta)
                call adapt_move(type2_spheres, iStep, type2_obs, type2_units%move_delta)                
            end if
            
            ! Observables writing
            write(type1_units%obsThermal, *) iStep, type1_obs%Epot, type1_obs%activ, &
                                             type1_obs%move_reject, type1_obs%rotate_reject
            write(type2_units%obsThermal, *) iStep, type2_obs%Epot, type2_obs%activ, &
                                             type2_obs%move_reject
            write(mix_obsThermal_unit, *) iStep, mix_Epot
            write(obsThermal_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot
            
            if (iStep == Nthermal) then ! Definite thermalised displacements
                write(output_unit, *) "Thermalisation : over"
                call type1_spheres%set_move_delta(type1_obs%move_rejectAvg, type1_units%report)
                call type1_spheres%set_rotate_delta(type1_obs%rotate_rejectAvg, type1_units%report)
                call type2_spheres%set_move_delta(type2_obs%move_rejectAvg, type2_units%report)
            end if
        
        else MC_Regime ! Thermalisation over -> Equilibrium
        
            ! Chemical potentials : Widom method
            call widom(type1_spheres, type1_obs, type2_spheres, mix)
            call widom(type2_spheres, type2_obs, type1_spheres, mix)
        
            ! Observables accumulations
            call type1_obs%accumulate()
            call type2_obs%accumulate()
            mix_EpotSum = mix_EpotSum + mix_Epot

            write(type1_units%obsEquilib, *) iStep, type1_obs%Epot, type1_obs%activ, &
                                             type1_obs%move_reject, type1_obs%rotate_reject
            write(type2_units%obsEquilib, *) iStep, type2_obs%Epot, type2_obs%activ, &
                                             type2_obs%move_reject
            write(mix_obsEquilib_unit, *) iStep, mix_Epot
            write(obsEquilib_unit, *) iStep, type1_obs%Epot + type2_obs%Epot + mix_Epot

            if (snap) then ! Snap shots of the configuration
                call type1_spheres%snap_positions(iStep, type1_units%snap_positions)
                call type1_spheres%snap_orientations(iStep, type1_units%snap_orientations)
                call type2_spheres%snap_positions(iStep, type2_units%snap_positions)
            end if
            
        end if MC_Regime
        
        ! Reinitialize quantities to prevent them from drifting.
        if (modulo(iStep, type1_spheres%get_reInit_iStep()) == 0) then
            call type1_spheres%Epot_reci_reInit_structure(iStep, type1_units%structure_modulus)
            call type1_spheres%reInit_totalMoment(iStep, type1_units%totalMoment_modulus)
        end if
    
    end do MC_Cycle
    call cpu_time(tFin)

    write(output_unit, *) "End of cycles"

! End ----------------------------------------------------------------------------------------------

    ! Tests & results

    call final(type1_spheres, type1_units, type1_obs)
    call final(type2_spheres, type2_units, type2_obs)
    call mix_final(mix, type1_spheres, type2_spheres, mix_report_unit, mix_Epot, mix_EpotSum, &
                   mix_Epot_conf)
    
    Epot = type1_obs%Epot + type2_obs%Epot + mix_Epot
    Epot_conf = type1_spheres%Epot_conf() + type2_spheres%Epot_conf() + mix_Epot_conf
    call test_consist(Epot, Epot_conf, report_unit)
    EpotSum = type1_obs%EpotSum + type2_obs%EpotSum + mix_EpotSum
    call print_results(Ncol, EpotSum, tFin-tIni, report_unit)
    
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
