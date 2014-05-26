!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program monte_carlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use module_types_micro, only: Monte_Carlo_Arguments
use module_monte_carlo_arguments, only: read_arguments
use class_physical_system

implicit none
    
    type(Physical_System) :: sys
    type(Monte_Carlo_Arguments) :: args
    type(json_file) :: json
    
    integer :: iStep

    call json_initialize()
    call json%load_file(filename = "data.json")
    call read_arguments(args)
    call sys%construct(json)
    call sys%init(json, args)
        
    write(output_unit, *) "Beginning of cycles"
    
    call sys%set_time_start()
    MC_Cycle: do iStep = 1, sys%get_num_thermalisation_steps() + sys%get_num_equilibrium_steps()
    
        call sys%random_changes()
        call sys%update_rejections()
        
        MC_Regime: if (iStep <= sys%get_num_thermalisation_steps()) then
            
            call sys%adapt_changes(iStep)
            call sys%write_observables_thermalisation(iStep)
            
            if (iStep == sys%get_num_thermalisation_steps()) then
                write(output_unit, *) "Thermalisation: over"
                call sys%fix_changes()
            end if
        
        else MC_Regime
        
            call sys%measure_chemical_potentials()
            call sys%accumulate_observables()
            call sys%write_observables_equilibrium(iStep)
            call sys%take_snapshots(iStep)
            
        end if MC_Regime
        
        call sys%reinitialize_quantites(iStep)
    
    end do MC_Cycle
    call sys%set_time_end()

    write(output_unit, *) "End of cycles"
    
    call sys%final(json)
    call sys%destroy()
    call json%destroy()
    
end program monte_carlo_canonical_bulk
