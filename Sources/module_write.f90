!> \brief Subroutines as tools for the main program

module module_write

use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_constants, only: PI, sigma3d
use data_box, only: Ndim, Lsize, Kmax
use data_monteCarlo, only: Temperature, Nstep, decorrelFactor, Nthermal
use data_potential, only: write_potential

implicit none

private
public open_units, mix_open_units, write_report, write_results, mix_write_results

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
    
    !> Total: write_report
    
    subroutine write_report(Ncol, Nmove, Nswitch, Nrotate, reset_iStep, report_unit)
    
        integer, intent(in) :: Ncol, Nmove, Nswitch, Nrotate, reset_iStep
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
        
        write(report_unit, *) "    reset_iStep = ", reset_iStep
        write(report_unit, *) "    write_potential = ", write_potential
    
    end subroutine write_report
    
    !> Total: Results
    
    subroutine write_results(Ncol, EpotSum, switch_rejectSum, duration, report_unit)
    
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
    
    end subroutine write_results
    
    !> Mix: Results
    
    subroutine mix_write_results(EpotSum, report_unit)
    
        real(DP), intent(in) :: EpotSum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", EpotSum/real(Nstep, DP)
    
    end subroutine mix_write_results
    
end module module_write
