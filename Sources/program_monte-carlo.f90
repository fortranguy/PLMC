!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program monteCarlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only : output_unit
use data_precisions, only : DP
use data_monteCarlo, only : decorrelFactor, Nthermal, Nadapt, Nstep
use data_distribution, only : snap
use class_neighbourCells
use class_hardSpheres
use class_interactingSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units
use module_algorithms, only : move, widom, rotate
use module_tools, only : init_randomSeed, set_initialCondition, print_report, test_consist, &
                         print_results, mix_print_results

implicit none
    
    ! Monte-Carlo variables
    integer :: iStep !< step counter
    integer :: iChange !< change counters
    integer :: iChangeRand !< random change
    integer :: Ncol !< Number of particles
    integer :: Nmove !< Number of displacements
    integer :: Nrotate !< Number of rotations
    integer :: iColRand !< random particle
    real(DP) :: random !< random number between 0 and 1
    real(DP) :: tIni, tFin !< CPU initial and final time
    
    ! Total physical system variables
    real(DP) :: Epot, EpotSum !< potential energy : at a Monte-Carlo step
    real(DP) :: Epot_conf !< potential energy : complete calculation from a configuration
    integer :: report_unit !< data & results file
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

    call type1_spheres%construct()
    call type2_spheres%construct()
    call mix%construct(type1_spheres%get_rMin(), type2_spheres%get_rMin())
    
    Ncol = type1_spheres%get_Ncol() + type2_spheres%get_Ncol()
    Nmove = decorrelFactor * Ncol
    Nrotate = decorrelFactor * type1_spheres%get_Ncol()
    
! Beginning ----------------------------------------------------------------------------------------
    
    write(output_unit, *) "Monte-Carlo Simulation : Canonical ensemble"

    ! Initialisations & reports
    
    open(newunit=report_unit, recl=4096, file="report.txt", status='new', action='write')
    open(newunit=obsThermal_unit, recl=4096, file="obsThermal.out", status='new', &
         action='write')
    open(newunit=obsEquilib_unit, recl=4096, file="obsEquilib.out", status='new', &
         action='write')
    write(obsEquilib_unit, *) "#", 1 ! 1 observable : energy
    call print_report(Ncol, Nmove, Nrotate, report_unit)
    !call init_randomSeed(report_unit)
    
    mix_EpotSum = 0._DP
    open(newunit=mix_report_unit, recl=4096, file="mix_report.txt", status='new', action='write')
    open(newunit=mix_Epot_unit, recl=4096, file="mix_Epot.tmp", status='new', action='write')
    open(newunit=mix_obsThermal_unit, recl=4096, file="mix_obsThermal.out", &
         status='new', action='write')
    open(newunit=mix_obsEquilib_unit, recl=4096, file="mix_obsEquilib.out", status='new', &
         action='write')
    write(mix_obsEquilib_unit, *) "#", 1 ! 1 observable : energy
    call mix%print_report(mix_report_unit)
    
    call type1_obs%init()
    call type1_units%open(type1_spheres%get_name())
    write(type1_units%obsEquilib, *) "#", 4 ! 4 observables
    call type1_spheres%print_density(Ncol, type1_units%report)
    
    call type2_obs%init()
    call type2_units%open(type2_spheres%get_name())
    write(type2_units%obsEquilib, *) "#", 3 ! 3 observables
    call type2_spheres%print_density(Ncol, type2_units%report)
    
    ! Initial condition
    
    call set_initialCondition(type1_spheres, type2_spheres, mix%get_rMin(), report_unit)
    
    call mix%test_overlap(type1_spheres, type2_spheres)
    call mix%Epot_init()
    call mix%Epot_print(mix_Epot_unit)
    call mix%set_cell_size()
    mix_Epot = mix%Epot_conf(type1_spheres, type2_spheres)
    
    call type1_spheres%test_overlap()
    call type1_spheres%Epot_init()
    call type1_spheres%Epot_real_print(type1_units%Epot)
    call type1_spheres%Epot_reci_count_waveVectors(type1_units%waveVectors)
    type1_obs%Epot = type1_spheres%Epot_conf()
    call type1_spheres%snap_positions_data(type1_units%snap_positions)
    call type1_spheres%snap_positions(0, type1_units%snapIni_positions)
    call type1_spheres%snap_orientations(0, type1_units%snapIni_orientations)
    call type1_spheres%construct_cells(type2_spheres, mix%get_cell_size(), mix%get_rCut())
    call type1_spheres%print_report(type1_units%report)
    
    call type2_spheres%test_overlap()
    call type2_spheres%Epot_init()
    call type2_spheres%Epot_print(type2_units%Epot)
    type2_obs%Epot = type2_spheres%Epot_conf()
    call type2_spheres%snap_positions_data(type2_units%snap_positions)
    call type2_spheres%snap_positions(0, type2_units%snapIni_positions)
    call type2_spheres%construct_cells(type1_spheres, mix%get_cell_size(), mix%get_rCut())
    call type2_spheres%print_report(type2_units%report)
    
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
            
            if (iChangeRand <= Nmove) then ! change = move
            
                ! Randomly choosing the type
                call random_number(random)
                iColRand = int(random*real(Ncol, DP)) + 1
                
                ! Moving a particle :
                if (iColRand <= type1_spheres%get_Ncol()) then
                    call move(type1_spheres, type1_obs, type2_spheres, mix, mix_Epot)
                    type1_obs%Nmove = type1_obs%Nmove + 1
                else
                    call move(type2_spheres, type2_obs, type1_spheres, mix, mix_Epot)
                    type2_obs%Nmove = type2_obs%Nmove + 1
                end if
                
            else ! change = rotate
     
                call rotate(type1_spheres, type1_obs)
                type1_obs%Nrotate = type1_obs%Nrotate + 1
                
            end if
            
        end do MC_Change
        
        ! Rejections rates updates
        call type1_obs%update_rejections()
        call type2_obs%update_rejections()
        
        MC_Regime : if (iStep <= Nthermal) then ! Thermalisation
            
            ! Displacements adaptation
            if (mod(iStep, Nadapt) /= 0) then ! Rejections accumulation
            
                type1_obs%move_rejectAdapt = type1_obs%move_rejectAdapt + type1_obs%move_reject
                type1_obs%rotate_rejectAdapt = type1_obs%rotate_rejectAdapt + type1_obs%rotate_reject
                type2_obs%move_rejectAdapt = type2_obs%move_rejectAdapt + type2_obs%move_reject
                
            else ! Average & adaptation
            
                type1_obs%move_rejectAvg = type1_obs%move_rejectAdapt / real(Nadapt-1, DP)
                type1_obs%move_rejectAdapt = 0._DP
                call type1_spheres%adapt_move_delta(type1_obs%move_rejectAvg)
                write(type1_units%move_delta, *) iStep, type1_spheres%get_move_delta(), &
                                                 type1_obs%move_rejectAvg
                
                
                type1_obs%rotate_rejectAvg = type1_obs%rotate_rejectAdapt / real(Nadapt-1, DP)
                type1_obs%rotate_rejectAdapt = 0._DP
                call type1_spheres%adapt_rotate_delta(type1_obs%rotate_rejectAvg)
                write(type1_units%rotate_delta, *) iStep, type1_spheres%get_rotate_delta(), &
                                                   type1_obs%rotate_rejectAvg
                
                type2_obs%move_rejectAvg = type2_obs%move_rejectAdapt / real(Nadapt-1, DP)
                type2_obs%move_rejectAdapt = 0._DP
                call type2_spheres%adapt_move_delta(type2_obs%move_rejectAvg)
                write(type2_units%move_delta, *) iStep, type2_spheres%get_move_delta(), &
                                                 type2_obs%move_rejectAvg
                
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
        
        ! Reinitialize to prevent it from drifting.
        if (modulo(iStep, type1_spheres%get_reInit_iStep()) == 0) then
            call type1_spheres%Epot_reci_reInit_structure(iStep, type1_units%structure_modulus)
            call type1_spheres%reInit_totalMoment(iStep, type1_units%totalMoment_modulus)
        end if
    
    end do MC_Cycle
    call cpu_time(tFin)

    write(output_unit, *) "End of cycles"

! End ----------------------------------------------------------------------------------------------

    ! Tests & results

    call type1_spheres%test_overlap()
    call type1_spheres%Epot_init()
    call type1_spheres%test_consist(type1_obs%Epot, type1_units%report)
    call type1_spheres%snap_positions(0, type1_units%snapFin_positions)
    call type1_spheres%snap_orientations(0, type1_units%snapFin_orientations)
    call type1_obs%print_results(type1_units%report)
    
    call type2_spheres%test_overlap()
    call type2_spheres%Epot_init()
    call type2_spheres%test_consist(type2_obs%Epot, type2_units%report)
    call type2_spheres%snap_positions(0, type2_units%snapFin_positions)
    call type2_obs%print_results(type2_units%report)
    
    call mix%test_overlap(type1_spheres, type2_spheres)
    call mix%Epot_init()
    mix_Epot_conf = mix%Epot_conf(type1_spheres, type2_spheres)
    call test_consist(mix_Epot, mix_Epot_conf, mix_report_unit)
    call mix_print_results(mix_EpotSum, mix_report_unit)
    
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
