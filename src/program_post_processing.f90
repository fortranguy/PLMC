!> \brief Monte Carlo simulation in canonical ensemble for a mixture

program monte_post_processing

use, intrinsic :: iso_fortran_env, only: output_unit
use module_types_micro, only: System_Arguments
use class_system, only: System_Post_Processing
use module_arguments_system, only: read_arguments_post

implicit none
    
    type(System_Post_Processing) :: sys
    type(System_Arguments) :: args
    
    integer :: i_step
    logical :: coordinates_set

    call read_arguments_post(args)
    call sys%construct(args)
    call sys%init(args)
        
    write(output_unit, *) "Beginning of cycles"
    
    call sys%set_time_start()
    MC_Cycle: do i_step = sys%get_num_thermalisation_steps() + 1, &
                          sys%get_num_thermalisation_steps() + sys%get_num_equilibrium_steps()
                          
            call sys%set_coordinates(i_step, coordinates_set)
            
            if (coordinates_set) then
                if (sys%get_first_set()) call sys%set_first()
                call sys%measure_chemical_potentials()
                call sys%accumulate_observables()
                !call sys%write_observables(i_step)
                call sys%reset_quantites(i_step)
            end if
    
    end do MC_Cycle
    call sys%set_time_end()

    write(output_unit, *) "End of cycles"
    
    call sys%final()
    call sys%destroy()
    
end program monte_post_processing
