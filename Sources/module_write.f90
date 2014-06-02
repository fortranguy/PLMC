!> \brief Subroutines as tools for the main program

module module_write

use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_constants, only: PI, sigma3d
use data_box, only: Ndim

implicit none

private
public open_units, write_data, mix_open_units, write_results, mix_write_results

contains

    !> Total: open units
    
    subroutine open_units(report_unit, observables_thermalisation_unit, observables_equilibrium_unit)
    
        integer, intent(out) :: report_unit, observables_thermalisation_unit, observables_equilibrium_unit
    
        open(newunit=report_unit, recl=4096, file="report.txt", status='new', action='write')
        open(newunit=observables_thermalisation_unit, recl=4096, file="observables_thermalisation.out", &
             status='new', action='write')
        open(newunit=observables_equilibrium_unit, recl=4096, file="observables_equilibrium.out", &
             status='new', action='write')
        write(observables_equilibrium_unit, *) "#", 1 ! 1 observable: energy
        
    end subroutine open_units

    !> Data: low level

    subroutine write_data(report_unit)
    
        integer, intent(in) :: report_unit

        write(report_unit, *) "Data micro: "

        write(report_unit ,*) "    Precision = ", DP
        write(report_unit ,*) "    Real zero = ", real_zero
        write(report_unit ,*) "    I/O tiny = ", io_tiny
        write(report_unit ,*) "    Energy consistency tiny = ", consist_tiny

        write(report_unit ,*) "    Pi = ", PI
        write(report_unit ,*) "    Sigma3d = ", sigma3d

    end subroutine write_data

    !> Mix: open units
    
    subroutine mix_open_units(mix_report_unit, mix_potential_unit, mix_observables_thermalisation_unit, &
                              mix_observables_equilibrium_unit)
                                  
        integer, intent(out) :: mix_report_unit, mix_potential_unit, mix_observables_thermalisation_unit, &
                                mix_observables_equilibrium_unit
    
        open(newunit=mix_report_unit, recl=4096, file="mix_report.txt", status='new', action='write')
        open(newunit=mix_potential_unit, recl=4096, file="mix_potential.tmp", status='new', action='write')
        open(newunit=mix_observables_thermalisation_unit, recl=4096, file="mix_observables_thermalisation.out", &
             status='new', action='write')
        open(newunit=mix_observables_equilibrium_unit, recl=4096, file="mix_observables_equilibrium.out", &
             status='new', action='write')
        write(mix_observables_equilibrium_unit, *) "#", 1 ! 1 observable: energy
        
     end subroutine mix_open_units
    
    !> Total: Results
    
    subroutine write_results(num_particles, num_equilibrium_steps, potential_sum, switch_rejectSum, duration, report_unit)
    
        integer, intent(in) :: num_particles, num_equilibrium_steps
        real(DP), intent(in) :: potential_sum, switch_rejectSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", potential_sum/real(num_equilibrium_steps, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   potential_sum/real(num_equilibrium_steps, DP)/real(num_particles, DP)
        write(report_unit, *) "    switch rejection rate = ", switch_rejectSum/real(num_equilibrium_steps, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine write_results
    
    !> Mix: Results
    
    subroutine mix_write_results(num_equilibrium_steps, potential_sum, report_unit)
    
        integer, intent(in) :: num_equilibrium_steps
        real(DP), intent(in) :: potential_sum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", potential_sum/real(num_equilibrium_steps, DP)
    
    end subroutine mix_write_results
    
end module module_write
