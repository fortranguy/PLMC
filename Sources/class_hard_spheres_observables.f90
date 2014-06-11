!> \brief Description of the Hard_Spheres_Observables class

module class_hard_spheres_observables

use data_precisions, only: DP

implicit none
private

    type, public :: Hard_Spheres_Observables
    
        ! Move
        integer :: move_num_hits = 0
        integer :: move_num_rejections = 0
        real(DP) :: move_rejection_rate = 0._DP
        real(DP) :: move_sum_rejection = 0._DP
        real(DP) :: move_rejection_adapt = 0._DP
        real(DP) :: move_rejection_average = 0._DP
        
        ! Potential energy
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
        
        ! Inverse of activity
        real(DP) :: inv_activity
        real(DP) :: sum_inv_activity = 0._DP
    
    contains
    
        procedure :: update_rejections => Hard_Spheres_Observables_update_rejections
        procedure :: write => Hard_Spheres_Observables_write
        procedure :: accumulate => Hard_Spheres_Observables_accumulate
        procedure :: write_results => Hard_Spheres_Observables_write_results
    
    end type Hard_Spheres_Observables
    
    type, extends(Hard_Spheres_Observables), public :: Dipolar_Hard_Spheres_Observables
        
        ! Rotate
        integer :: rotate_num_hits = 0
        integer :: rotate_num_rejections = 0
        real(DP) :: rotate_rejection_rate = 0._DP
        real(DP) :: rotate_sum_rejection = 0._DP
        real(DP) :: rotate_rejection_adapt = 0._DP
        real(DP) :: rotate_rejection_average = 0._DP
        
    end type Dipolar_Hard_Spheres_Observables
    
    type, public :: Between_Hard_Spheres_Observables
    
        real(DP) :: potential_energy
        real(DP) :: potential_energy_sum = 0._DP
        
    contains
    
        procedure :: write => Between_Hard_Spheres_Observables_write
        procedure :: accumulate => Between_Hard_Spheres_Observables_accumulate
    
    end type Between_Hard_Spheres_Observables
    
contains
    
    pure subroutine Hard_Spheres_Observables_update_rejections(this)
    
        class(Hard_Spheres_Observables), intent(inout) :: this
        
        this%move_rejection_rate = real(this%move_num_rejections, DP)/real(this%move_num_hits, DP)
        this%move_num_rejections = 0
        this%move_num_hits = 0
        
         select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                this%rotate_rejection_rate = real(this%rotate_num_rejections, DP) / &
                                             real(this%rotate_num_hits, DP)
                this%rotate_num_rejections = 0
                this%rotate_num_hits = 0
        end select
    
    end subroutine Hard_Spheres_Observables_update_rejections
    
    ! Write
    
    subroutine Hard_Spheres_Observables_write(this, i_step, observables_unit)
    
        class(Hard_Spheres_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(observables_unit, *) i_step, this%potential_energy, this%inv_activity, &
                                           this%move_rejection_rate, this%rotate_rejection_rate
            class default
                write(observables_unit, *) i_step, this%potential_energy, this%inv_activity, &
                                           this%move_rejection_rate
        end select
    
    end subroutine Hard_Spheres_Observables_write
    
    subroutine Between_Hard_Spheres_Observables_write(this, i_step, observables_unit)
    
        class(Between_Hard_Spheres_Observables), intent(in) :: this
        integer, intent(in) :: i_step, observables_unit
        
        write(observables_unit, *) i_step, this%potential_energy
    
    end subroutine Between_Hard_Spheres_Observables_write
    
    !> Accumulations
    
    pure subroutine Hard_Spheres_Observables_accumulate(this)
    
        class(Hard_Spheres_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
        this%sum_inv_activity = this%sum_inv_activity + this%inv_activity
        this%move_sum_rejection = this%move_sum_rejection + this%move_rejection_rate
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                this%rotate_sum_rejection = this%rotate_sum_rejection + this%rotate_rejection_rate
        end select
    
    end subroutine Hard_Spheres_Observables_accumulate
    
    pure subroutine Between_Hard_Spheres_Observables_accumulate(this)
    
        class(Between_Hard_Spheres_Observables), intent(inout) :: this
        
        this%potential_energy_sum = this%potential_energy_sum + this%potential_energy
    
    end subroutine Between_Hard_Spheres_Observables_accumulate
    
    !> Results
    
    subroutine Hard_Spheres_Observables_write_results(this, temperature, num_equilibrium_steps, &
                                                      report_unit)

        class(Hard_Spheres_Observables), intent(in) :: this
        real(DP), intent(in) :: temperature
        integer, intent(in) :: num_equilibrium_steps
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiEx
            
        write(report_unit, *) "Results: "
        
        write(report_unit, *) "    average energy = ", this%potential_energy_sum / &
                                                       real(num_equilibrium_steps, DP)
        potChiEx = -temperature*log(this%sum_inv_activity/real(num_equilibrium_steps, DP))
        write(report_unit, *) "    average excess chemical potential = ", potChiEx
        write(report_unit, *) "    move rejection rate = ", &
                                   this%move_sum_rejection/real(num_equilibrium_steps, DP)
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(report_unit, *) "    rotate rejection rate = ", &
                                      this%rotate_sum_rejection/real(num_equilibrium_steps, DP)
        end select
    
    end subroutine Hard_Spheres_Observables_write_results

end module class_hard_spheres_observables
