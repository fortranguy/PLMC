!> \brief Subroutines as tools for the main program

module module_write

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_value, json_value_create, to_object, json_value_add

implicit none

private
public write_results

contains
    
    !> Total: Results
    
    subroutine write_results(num_particles, num_equilibrium_steps, potential_energy_sum, switch_rejectSum, &
                             duration, report_json)
    
        integer, intent(in) :: num_particles, num_equilibrium_steps
        real(DP), intent(in) :: potential_energy_sum, switch_rejectSum
        real(DP), intent(in) :: duration
        type(json_value), pointer, intent(in) :: report_json

        type(json_value), pointer :: results_json

        call json_value_create(results_json)
        call to_object(results_json, "Results")
        call json_value_add(report_json, results_json)

        call json_value_add(results_json, "average energy", &
                                          potential_energy_sum/real(num_equilibrium_steps, DP))
        call json_value_add(results_json, "average energy per particule", &
                                          potential_energy_sum / real(num_equilibrium_steps, DP) / &
                                          real(num_particles, DP)) ! DHS only?
        call json_value_add(results_json, "switch rejection rate", switch_rejectSum / &
                                          real(num_equilibrium_steps, DP))
        call json_value_add(results_json, "duration (min)", duration/60._DP)

        nullify(results_json)
    
    end subroutine write_results
    
end module module_write
