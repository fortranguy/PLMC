!> \brief Monte-Carlo simulation in canonical ensemble for a mixture

program monteCarlo_canonical_bulk

use, intrinsic :: iso_fortran_env, only: output_unit
use module_types, only: argument_seed, argument_initial
use module_monteCarlo_arguments, only: read_arguments
use class_physicalSystem

implicit none

    integer :: iStep
    
    type(PhysicalSystem) :: sys
    type(argument_seed) :: arg_seed
    type(argument_initial) :: arg_init

    call read_arguments(arg_seed, arg_init)
    call sys%construct()
    call sys%init(arg_seed, arg_init)
        
    write(output_unit, *) "Beginning of cycles"
    
    !call cpu_time(tIni)
    MC_Cycle: do iStep = 1, sys%get_Nthermal() + sys%get_Nstep()
    
        call sys%random_changes() 
        call sys%update_rejections()
        
        MC_Regime: if (iStep <= sys%get_Nthermal()) then ! Thermalisation
            
            call sys%adapt_changes(iStep)            
            call sys%write_observables_thermalisation(iStep)
            
            if (iStep == sys%get_Nthermal()) then
                write(output_unit, *) "Thermalisation: over"
                call sys%fix_changes()
            end if
        
        else MC_Regime ! Thermalisation over -> Equilibrium
        
            call sys%measure_chemical_potentials()        
            call sys%accumulate_observables()            
            call sys%write_observables_equilibrium(iStep)
            call sys%take_snapshots(iStep)
            
        end if MC_Regime
        
        call sys%reinitialize_quantites(iStep)
    
    end do MC_Cycle
    !call cpu_time(tFin)

    write(output_unit, *) "End of cycles"
    
    call sys%final()    
    call sys%destroy()
    
end program monteCarlo_canonical_bulk
