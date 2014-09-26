!> \brief Description of the Hard_Spheres_Monte_Carlo_Observables class

module class_hard_spheres_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_value, json_value_create, to_object, json_value_add
use class_discrete_observable, only: Discrete_Observables, Adapting_Discrete_Observables

implicit none
private

    type, public :: Hard_Spheres_Monte_Carlo_Observables
    
        type(Adapting_Discrete_Observables) :: move
        
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
    
    contains
    
        procedure :: write => Hard_Spheres_Monte_Carlo_Observables_write
        procedure :: accumulate => Hard_Spheres_Monte_Carlo_Observables_accumulate
        procedure :: write_results => Hard_Spheres_Monte_Carlo_Observables_write_results
    
    end type Hard_Spheres_Monte_Carlo_Observables
    
    type, extends(Hard_Spheres_Monte_Carlo_Observables), public :: Dipolar_Hard_Spheres_Monte_Carlo_Observables
        
        type(Adapting_Discrete_Observables) :: rotation
        
    end type Dipolar_Hard_Spheres_Monte_Carlo_Observables
    
    type, public :: Between_Hard_Spheres_Monte_Carlo_Observables
    
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
        
    contains
    
        procedure :: write => Between_Hard_Spheres_Monte_Carlo_Observables_write
        procedure :: accumulate => Between_Hard_Spheres_Monte_Carlo_Observables_accumulate
        procedure :: write_results => Between_Hard_Spheres_Monte_Carlo_Observables_write_results
    
    end type Between_Hard_Spheres_Monte_Carlo_Observables

    type, public :: Hard_Spheres_Post_Processing_Observables

        real(DP) :: inv_activity
        real(DP) :: sum_inv_activity = 0._DP

    contains

        procedure :: write => Hard_Spheres_Post_Processing_Observables_write
        procedure :: accumulate => Hard_Spheres_Post_Processing_Observables_accumulate
        procedure :: write_results => Hard_Spheres_Post_Processing_Observables_write_results
        
    end type Hard_Spheres_Post_Processing_Observables
    
    type, extends(Hard_Spheres_Post_Processing_Observables), public :: &
        Dipolar_Hard_Spheres_Post_Processing_Observables
        type(Discrete_Observables) :: local_field
    end type Dipolar_Hard_Spheres_Post_Processing_Observables
    
contains
    
    subroutine Hard_Spheres_Monte_Carlo_Observables_write(this, i_step, observables_unit)
    
        class(Hard_Spheres_Monte_Carlo_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Monte_Carlo_Observables)
                write(observables_unit, *) i_step, this%potential_energy, &
                                           this%move%rejection_rate, this%rotation%rejection_rate
            class default
                write(observables_unit, *) i_step, this%potential_energy, &
                                           this%move%rejection_rate
        end select
    
    end subroutine Hard_Spheres_Monte_Carlo_Observables_write
    
    subroutine Between_Hard_Spheres_Monte_Carlo_Observables_write(this, i_step, observables_unit)
    
        class(Between_Hard_Spheres_Monte_Carlo_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        write(observables_unit, *) i_step, this%potential_energy
    
    end subroutine Between_Hard_Spheres_Monte_Carlo_Observables_write

    subroutine Hard_Spheres_Post_Processing_Observables_write(this, i_step, observables_unit)

        class(Hard_Spheres_Post_Processing_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        select type (this)
            type is (Hard_Spheres_Post_Processing_Observables)
                write(observables_unit, *) i_step, this%inv_activity
            type is (Dipolar_Hard_Spheres_Post_Processing_Observables)
                write(observables_unit, *) i_step, this%inv_activity, this%local_field%rejection_rate
        end select

    end subroutine Hard_Spheres_Post_Processing_Observables_write
    
    pure subroutine Hard_Spheres_Monte_Carlo_Observables_accumulate(this)
    
        class(Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
        this%move%sum_rejection = this%move%sum_rejection + this%move%rejection_rate
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Monte_Carlo_Observables)
                this%rotation%sum_rejection = this%rotation%sum_rejection + &
                                              this%rotation%rejection_rate
        end select
    
    end subroutine Hard_Spheres_Monte_Carlo_Observables_accumulate
    
    pure subroutine Between_Hard_Spheres_Monte_Carlo_Observables_accumulate(this)
    
        class(Between_Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
    
    end subroutine Between_Hard_Spheres_Monte_Carlo_Observables_accumulate

    pure subroutine Hard_Spheres_Post_Processing_Observables_accumulate(this)

        class(Hard_Spheres_Post_Processing_Observables), intent(inout) :: this

        this%sum_inv_activity = this%sum_inv_activity + this%inv_activity
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Post_Processing_Observables)            
                this%local_field%sum_rejection = this%local_field%sum_rejection + &
                                                 this%local_field%rejection_rate
        end select

    end subroutine Hard_Spheres_Post_Processing_Observables_accumulate
    
    subroutine Hard_Spheres_Monte_Carlo_Observables_write_results(this, num_equilibrium_steps, &
                                                                  report_json)

        class(Hard_Spheres_Monte_Carlo_Observables), intent(in) :: this
        integer, intent(in) :: num_equilibrium_steps
        type(json_value), pointer, intent(in) :: report_json
        
        type(json_value), pointer :: results_json

        call json_value_create(results_json)
        call to_object(results_json, "Results")
        call json_value_add(report_json, results_json)

        call json_value_add(results_json, "average energy", &
                                          this%potential_energy_sum / real(num_equilibrium_steps, DP))
        call json_value_add(results_json, "move rejection rate", &
                                          this%move%sum_rejection / real(num_equilibrium_steps, DP))

        select type (this)
            type is (Dipolar_Hard_Spheres_Monte_Carlo_Observables)
                call json_value_add(results_json, "rotation rejection rate", &
                                                   this%rotation%sum_rejection / &
                                                   real(num_equilibrium_steps, DP))
        end select

        nullify(results_json)
    
    end subroutine Hard_Spheres_Monte_Carlo_Observables_write_results
    
    subroutine Between_Hard_Spheres_Monte_Carlo_Observables_write_results(this, num_equilibrium_steps, &
                                                                          report_json)
                                          
        class(Between_Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: this
        integer, intent(in) :: num_equilibrium_steps
        type(json_value), pointer, intent(in) :: report_json

        type(json_value), pointer :: results_json

        call json_value_create(results_json)
        call to_object(results_json, "Results")
        call json_value_add(report_json, results_json)

        call json_value_add(results_json, "average energy", &
                                          this%potential_energy_sum / real(num_equilibrium_steps, DP))
        
        nullify(results_json)
    
    end subroutine Between_Hard_Spheres_Monte_Carlo_Observables_write_results

    subroutine Hard_Spheres_Post_Processing_Observables_write_results(this, temperature, &
                                                                      num_equilibrium_steps, &
                                                                      widom_num_particles, &
                                                                      report_json)

        class(Hard_Spheres_Post_Processing_Observables), intent(in) :: this
        real(DP), intent(in) :: temperature
        integer, intent(in) :: num_equilibrium_steps
        integer, intent(in) :: widom_num_particles
        type(json_value), pointer, intent(in) :: report_json

        real(DP) :: chemical_potential_excess
        type(json_value), pointer :: results_json

        call json_value_create(results_json)
        call to_object(results_json, "Results")
        call json_value_add(report_json, results_json)

        if (widom_num_particles > 0) then
            chemical_potential_excess = -temperature * log(this%sum_inv_activity / &
                                                    real(num_equilibrium_steps, DP))
            call json_value_add(results_json, "average excess chemical potential", &
                                              chemical_potential_excess)
        end if
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Post_Processing_Observables)
                call json_value_add(results_json, "local field rejection rate", &
                                    this%local_field%sum_rejection / real(num_equilibrium_steps, DP))
        end select

        nullify(results_json)

    end subroutine Hard_Spheres_Post_Processing_Observables_write_results

end module class_hard_spheres_observables
