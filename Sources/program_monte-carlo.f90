!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program monteCarlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only: output_unit
use module_types, only: monteCarlo_arguments
use module_monteCarlo_arguments, only: read_arguments
use class_physicalSystem

implicit none
    
    type(PhysicalSystem) :: sys
    type(monteCarlo_arguments) :: args
    
    integer :: iStep

    call read_arguments(args)
    call sys%construct()
    call sys%init(args)
        
    write(output_unit, *) "Beginning of cycles"
    
    call sys%set_time_start()
    MC_Cycle: do iStep = 1, sys%get_Nthermal() + sys%get_Nstep()
    
        call sys%random_changes()
        call sys%update_rejections()
        
        MC_Regime: if (iStep <= sys%get_Nthermal()) then
            
            call sys%adapt_changes(iStep)
            call sys%write_observables_thermalisation(iStep)
            
            if (iStep == sys%get_Nthermal()) then
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
    
    call sys%final()
    call sys%destroy()
    
end program monteCarlo_canonical_bulk
