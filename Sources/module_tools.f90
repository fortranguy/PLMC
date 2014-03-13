!> \brief Subroutines as tools for the main program

module module_tools

use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_constants, only: PI, sigma3d
use data_box, only: Ndim, Lsize, Kmax
use data_monteCarlo, only: Temperature, Nadapt, Nstep, decorrelFactor, Nthermal
use data_potential, only: print_potential
use module_types, only: argument_seed, argument_initial
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables
use class_units

implicit none

private
public open_units, mix_open_units, &
       print_report, mix_init, mix_final, init, final, adapt_move, adapt_rotate, test_consist, &
       print_results, mix_print_results

contains

    !> Total: open units
    
    subroutine open_units(report_unit, obsThermal_unit, obsEquilib_unit)
    
        integer, intent(out) :: report_unit, obsThermal_unit, obsEquilib_unit
    
        open(newunit=report_unit, recl=4096, file="report.txt", status='new', action='write')
        open(newunit=obsThermal_unit, recl=4096, file="obsThermal.out", status='new', &
             action='write')
        open(newunit=obsEquilib_unit, recl=4096, file="obsEquilib.out", status='new', &
             action='write')
        write(obsEquilib_unit, *) "#", 1 ! 1 observable: energy
        
    end subroutine open_units

    !> Mix: open units
    
    subroutine mix_open_units(mix_report_unit, mix_Epot_unit, mix_obsThermal_unit, &
                                  mix_obsEquilib_unit)
                                  
        integer, intent(out) :: mix_report_unit, mix_Epot_unit, mix_obsThermal_unit, &
                                  mix_obsEquilib_unit
    
        open(newunit=mix_report_unit, recl=4096, file="mix_report.txt", status='new', action='write')
        open(newunit=mix_Epot_unit, recl=4096, file="mix_Epot.tmp", status='new', action='write')
        open(newunit=mix_obsThermal_unit, recl=4096, file="mix_obsThermal.out", &
             status='new', action='write')
        open(newunit=mix_obsEquilib_unit, recl=4096, file="mix_obsEquilib.out", status='new', &
             action='write')
        write(mix_obsEquilib_unit, *) "#", 1 ! 1 observable: energy
        
     end subroutine mix_open_units
    
    !> Total: print_report
    
    subroutine print_report(Ncol, Nmove, Nswitch, Nrotate, report_unit)
    
        integer, intent(in) :: Ncol, Nmove, Nswitch, Nrotate
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data: "
        
        write(report_unit ,*) "    Precision = ", DP
        write(report_unit ,*) "    Real zero = ", real_zero
        write(report_unit ,*) "    I/O tiny = ", io_tiny
        write(report_unit ,*) "    Energy consistency tiny = ", consist_tiny
        
        write(report_unit ,*) "    Pi = ", PI
        write(report_unit ,*) "    Sigma3d = ", sigma3d
        
        write(report_unit ,*) "    Lsize(:) = ", Lsize(:)
        write(report_unit ,*) "    Volume = ", product(Lsize)
        write(report_unit ,*) "    Kmax(:) = ", Kmax(:)
        write(report_unit ,*) "    NwaveVectors =", (2*Kmax(1)+1) * (2*Kmax(2)+1) * (2*Kmax(3)+1)
        write(report_unit ,*) "    Ncol = ", Ncol
        write(report_unit ,*) "    Temperature = ", Temperature
        
        write(report_unit, *) "    Nstep = ", Nstep
        write(report_unit, *) "    Nthermal = ", Nthermal
        write(report_unit, *) "    decorrelFactor = ", decorrelFactor
        write(report_unit, *) "    Nmove = ", Nmove
        write(report_unit, *) "    Nswitch = ", Nswitch
        write(report_unit, *) "    Nrotate = ", Nrotate
    
    end subroutine print_report
    
    !> Mix initialisation
    
    subroutine mix_init(mix, type1, type2, mix_Epot_unit, mix_Epot)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_Epot_unit
        real(DP), intent(out) :: mix_Epot
    
        call mix%test_overlap(type1, type2)
        call mix%Epot_init()
        if (print_potential) then
            call mix%Epot_print(mix_Epot_unit)
        end if
        call mix%set_cell_size()
        mix_Epot = mix%Epot_conf(type1, type2)
    
    end subroutine mix_init
    
    !> Mix finalization
    
    subroutine mix_final(mix, type1, type2, mix_report_unit, mix_Epot, mix_EpotSum, mix_Epot_conf)
    
        class(MixingPotential), intent(inout) :: mix
        class(HardSpheres), intent(in) :: type1, type2
        integer, intent(in) :: mix_report_unit
        real(DP), intent(in) :: mix_Epot, mix_EpotSum
        real(DP), intent(out) :: mix_Epot_conf
        
        call mix%test_overlap(type1, type2)
        call mix%Epot_init()
        mix_Epot_conf = mix%Epot_conf(type1, type2)
        call test_consist(mix_Epot, mix_Epot_conf, mix_report_unit)
        call mix_print_results(mix_EpotSum, mix_report_unit)
    
    end subroutine mix_final
    
    !> Spheres initialisations
    
    subroutine init(this, other, mix, this_units, this_Epot)
    
        class(HardSpheres), intent(inout) :: this
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        class(Units), intent(in) :: this_units
        real(DP), intent(inout) :: this_Epot
        
        call this%test_overlap()
        call this%snap_data(this_units%snap_positions)
        call this%snap_positions(0, this_units%snapIni_positions)
        call this%Epot_init()
        
        if (print_potential) then
            call this%Epot_print(this_units%Epot)
        end if
        select type (this)
            type is (DipolarSpheres)
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_data(this_units%snap_orientations)
                        call this%snap_orientations(0, this_units%snapIni_orientations)
                        if (print_potential) then
                            call this%Epot_real_print(this_units%Epot_real)
                        end if
                        call this%Epot_reci_count_waveVectors(this_units%waveVectors)
                end select
        end select
        this_Epot = this%Epot_conf()
        
        call this%construct_cells(other, mix%get_cell_size(), mix%get_rCut())
        call this%print_report(this_units%report)
    
    end subroutine init
    
    !> Spheres finalizations
    
    subroutine final(this, this_units, this_obs)
    
        class(HardSpheres), intent(inout) :: this
        class(Units), intent(in) :: this_units
        class(Observables), intent(in) :: this_obs
        
        call this%test_overlap()
        call this%Epot_init()
        call test_consist(this_obs%Epot, this%Epot_conf(), this_units%report)
        call this%snap_positions(0, this_units%snapFin_positions)
        call this_obs%print_results(this_units%report)
        
        select type (this)
            type is (DipolarSpheres)
                select type (this_units)
                    type is (MoreUnits)
                        call this%snap_orientations(0, this_units%snapFin_orientations)
                end select
        end select
    
    end subroutine final
    
    !> Change: average & adaptation
    
    subroutine adapt_move(this, iStep, obs, move_unit)
    
        class(HardSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(Observables), intent(inout) :: obs
        integer, intent(in) :: move_unit
    
        obs%move_rejectAvg = obs%move_rejectAdapt / real(Nadapt-1, DP)
        obs%move_rejectAdapt = 0._DP
        call this%adapt_move_delta(obs%move_rejectAvg)
        write(move_unit, *) iStep, this%get_move_delta(), obs%move_rejectAvg
    
    end subroutine adapt_move
    
    subroutine adapt_rotate(this, iStep, obs, rotate_unit)
    
        class(DipolarSpheres), intent(inout) :: this
        integer, intent(in) :: iStep
        class(MoreObservables), intent(inout) :: obs
        integer, intent(in) :: rotate_unit
        
        obs%rotate_rejectAvg = obs%rotate_rejectAdapt / real(Nadapt-1, DP)
        obs%rotate_rejectAdapt = 0._DP
        call this%adapt_rotate_delta(obs%rotate_rejectAvg)
        write(rotate_unit, *) iStep, this%get_rotate_delta(), obs%rotate_rejectAvg
        
    end subroutine adapt_rotate
    
    !> Total & Mix: consistency test
    
    subroutine test_consist(Epot, Epot_conf, report_unit)
    
        real(DP), intent(in) :: Epot, Epot_conf
        integer, intent(in) :: report_unit
        
        real(DP) :: difference
        
        write(report_unit, *) "Consistency test: "
        write(report_unit, *) "    Epot = ", Epot
        write(report_unit, *) "    Epot_conf = ", Epot_conf
        
        if (abs(Epot_conf) < real_zero) then
            difference = abs(Epot_conf-Epot)
            write(report_unit, *) "    absolute difference = ", difference
        else
            difference = abs((Epot_conf-Epot)/Epot_conf)
            write(report_unit, *) "    relative difference = ", difference
        end if
        
        if (difference > consist_tiny) then ! not sufficient for HS ?
            write(report_unit, *) "    WARNING !"
        else
            write(report_unit, *) "    OK !"
        end if
    
    end subroutine test_consist
    
    !> Total: Results
    
    subroutine print_results(Ncol, EpotSum, switch_rejectSum, duration, report_unit)
    
        integer, intent(in) :: Ncol
        real(DP), intent(in) :: EpotSum, switch_rejectSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   EpotSum/real(Nstep, DP)/real(Ncol, DP)
        write(report_unit, *) "    switch rejection rate = ", switch_rejectSum/real(Nstep, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine print_results
    
    !> Mix: Results
    
    subroutine mix_print_results(EpotSum, report_unit)
    
        real(DP), intent(in) :: EpotSum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
    
    end subroutine mix_print_results
    
end module module_tools
