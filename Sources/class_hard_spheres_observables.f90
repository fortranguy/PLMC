!> \brief Description of the Hard_Spheres_Observables class

module class_hard_spheres_observables

use data_precisions, only: DP
use class_discrete_observable, only: Adapting_Discrete_Observables

implicit none
private

    type, public :: Hard_Spheres_Observables
    
        type(Adapting_Discrete_Observables) :: move
        
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
        
        real(DP) :: inv_activity
        real(DP) :: sum_inv_activity = 0._DP
    
    contains
    
        procedure :: write => Hard_Spheres_Observables_write
        procedure :: accumulate => Hard_Spheres_Observables_accumulate
        procedure :: write_results => Hard_Spheres_Observables_write_results
    
    end type Hard_Spheres_Observables
    
    type, extends(Hard_Spheres_Observables), public :: Dipolar_Hard_Spheres_Observables
        
        type(Adapting_Discrete_Observables) :: rotation
        
    end type Dipolar_Hard_Spheres_Observables
    
    type, public :: Between_Hard_Spheres_Observables
    
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
        
    contains
    
        procedure :: write => Between_Hard_Spheres_Observables_write
        procedure :: accumulate => Between_Hard_Spheres_Observables_accumulate
    
    end type Between_Hard_Spheres_Observables
    
contains
    
    subroutine Hard_Spheres_Observables_write(this, i_step, observables_unit)
    
        class(Hard_Spheres_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(observables_unit, *) i_step, this%potential_energy, this%inv_activity, &
                                           this%move%rejection_rate, this%rotation%rejection_rate
            class default
                write(observables_unit, *) i_step, this%potential_energy, this%inv_activity, &
                                           this%move%rejection_rate
        end select
    
    end subroutine Hard_Spheres_Observables_write
    
    subroutine Between_Hard_Spheres_Observables_write(this, i_step, observables_unit)
    
        class(Between_Hard_Spheres_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        write(observables_unit, *) i_step, this%potential_energy
    
    end subroutine Between_Hard_Spheres_Observables_write
    
    pure subroutine Hard_Spheres_Observables_accumulate(this)
    
        class(Hard_Spheres_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
        this%sum_inv_activity = this%sum_inv_activity + this%inv_activity
        this%move%sum_rejection = this%move%sum_rejection + this%move%rejection_rate
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                this%rotation%sum_rejection = this%rotation%sum_rejection + &
                                              this%rotation%rejection_rate
        end select
    
    end subroutine Hard_Spheres_Observables_accumulate
    
    pure subroutine Between_Hard_Spheres_Observables_accumulate(this)
    
        class(Between_Hard_Spheres_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
    
    end subroutine Between_Hard_Spheres_Observables_accumulate
    
    subroutine Hard_Spheres_Observables_write_results(this, temperature, num_equilibrium_steps, &
                                                      report_unit)

        class(Hard_Spheres_Observables), intent(in) :: this
        real(DP), intent(in) :: temperature
        integer, intent(in) :: num_equilibrium_steps
        integer, intent(in) :: report_unit
        
        real(DP) :: chemical_potential_excess
            
        write(report_unit, *) "Results: "
        
        write(report_unit, *) "    average energy = ", this%potential_energy_sum / &
                                                       real(num_equilibrium_steps, DP)
        chemical_potential_excess = -temperature * log(this%sum_inv_activity / &
                                    real(num_equilibrium_steps, DP))
        write(report_unit, *) "    average excess chemical potential = ", chemical_potential_excess
        write(report_unit, *) "    move rejection rate = ", &
                                   this%move%sum_rejection/real(num_equilibrium_steps, DP)
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(report_unit, *) "    rotation rejection rate = ", &
                                      this%rotation%sum_rejection/real(num_equilibrium_steps, DP)
        end select
    
    end subroutine Hard_Spheres_Observables_write_results

end module class_hard_spheres_observables
