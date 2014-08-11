!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program monte_carlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use json_module, only: json_file, json_initialize
use module_types_micro, only: Monte_Carlo_Arguments
use class_physical_system, only: Physical_System
use module_arguments_monte_carlo, only: read_arguments

implicit none
    
    type(Physical_System) :: sys
    type(Monte_Carlo_Arguments) :: args
    type(json_file) :: data_json
    character(len=*), parameter :: data_filename = "data.json"
    logical :: data_exists
    
    integer :: i_step

    inquire(file=data_filename, exist=data_exists)
    if (.not. data_exists) then
        write(error_unit, *) data_filename, " doesn't exist."
        error stop
    end if

    call json_initialize()
    call data_json%load_file(filename = data_filename)
    call read_arguments(args)
    call sys%construct(data_json, args)
    call sys%init(data_json, args)
        
    write(output_unit, *) "Beginning of cycles"
    
    call sys%set_time_start()
    MC_Cycle: do i_step = 1, sys%get_num_thermalisation_steps() + sys%get_num_equilibrium_steps()
    
        call sys%random_changes()
        call sys%update_rejections()
        
        MC_Regime: if (i_step <= sys%get_num_thermalisation_steps()) then
            
            call sys%adapt_changes(i_step)
            call sys%write_observables_thermalisation(i_step)
            
            if (i_step == sys%get_num_thermalisation_steps()) then
                write(output_unit, *) "Thermalisation: over"
                call sys%fix_changes()
            end if
        
        else MC_Regime
        
            call sys%measure_chemical_potentials()
            call sys%accumulate_observables()
            call sys%write_observables_equilibrium(i_step)
            call sys%take_snapshots(i_step)
            
        end if MC_Regime
        
        call sys%reinitialize_quantites(i_step)
    
    end do MC_Cycle
    call sys%set_time_end()

    write(output_unit, *) "End of cycles"
    
    call sys%final(data_json)
    call sys%destroy()
    call data_json%destroy()
    
end program monte_carlo_canonical_bulk
