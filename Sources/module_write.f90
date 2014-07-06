!> \brief Subroutines as tools for the main program

module module_write

use data_precisions, only: DP, real_zero, io_tiny, consist_tiny
use data_constants, only: PI, sigma3d
use data_box, only: num_dimensions
use module_types_micro, only: Box_Parameters
use class_hard_spheres, only: Hard_Spheres

implicit none

private
public open_units, write_data, write_spheres_density, write_results, between_spheres_write_results

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
        
        write(report_unit ,*) "    Number of dimensions = ", num_dimensions

        write(report_unit ,*) "    Pi = ", PI
        write(report_unit ,*) "    Sigma3d = ", sigma3d

    end subroutine write_data
    
    !> Total: Results
    
    subroutine write_results(num_particles, num_equilibrium_steps, potential_energy_sum, switch_rejectSum, &
                             duration, report_unit)
    
        integer, intent(in) :: num_particles, num_equilibrium_steps
        real(DP), intent(in) :: potential_energy_sum, switch_rejectSum
        real(DP), intent(in) :: duration
        integer, intent(in) :: report_unit
            
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", potential_energy_sum/real(num_equilibrium_steps, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   potential_energy_sum / real(num_equilibrium_steps, DP) / &
                                   real(num_particles, DP)
        write(report_unit, *) "    switch rejection rate = ", switch_rejectSum / &
                                   real(num_equilibrium_steps, DP)
        write(report_unit, *) "    duration =", duration/60._DP, "min"
    
    end subroutine write_results
    
    !> Mix: Results
    
    subroutine between_spheres_write_results(num_equilibrium_steps, potential_energy_sum, report_unit)
    
        integer, intent(in) :: num_equilibrium_steps
        real(DP), intent(in) :: potential_energy_sum
        integer, intent(in) :: report_unit
    
        write(report_unit, *) "Results: "
        write(report_unit, *) "    average energy = ", potential_energy_sum/real(num_equilibrium_steps, DP)
    
    end subroutine between_spheres_write_results
    
    !> Write density and compacity
    
    subroutine write_spheres_density(Box, this_spheres, report_unit)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres
        integer, intent(in) :: report_unit
        
        real(DP) :: density, compacity, concentration
        
        density = real(this_spheres%get_num_particles() + 1, DP) / product(Box%size) ! cheating ? cf. Widom
        compacity = 4._DP/3._DP*PI*(this_spheres%get_diameter()/2._DP)**3 * density
        concentration = real(this_spheres%get_num_particles(), DP) / real(Box%num_particles, DP)
        
        write(report_unit, *) "Density: "
        write(report_unit, *) "    density = ", density
        write(report_unit, *) "    compacity = ", compacity
        write(report_unit, *) "    concentration = ", concentration
    
    end subroutine write_spheres_density
    
end module module_write
