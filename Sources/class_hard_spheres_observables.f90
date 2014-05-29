!> \brief Description of the Hard_Spheres_Observables class

module class_hard_spheres_observables

use data_precisions, only: DP

implicit none
private

    type, public :: Hard_Spheres_Observables
    
        ! Move
        integer :: move_num_hits = 0
        integer :: move_Nreject = 0
        real(DP) :: move_reject = 0._DP
        real(DP) :: move_rejectSum = 0._DP
        real(DP) :: move_rejectAdapt = 0._DP
        real(DP) :: move_rejectAvg = 0._DP
        
        ! Potential energy
        real(DP) :: potential
        real(DP) :: potential_sum = 0._DP
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum = 0._DP
    
    contains
    
        procedure :: update_rejections => Hard_Spheres_Observables_update_rejections
        procedure :: write => Hard_Spheres_Observables_write
        procedure :: accumulate => Hard_Spheres_Observables_accumulate
        procedure :: write_results => Hard_Spheres_Observables_write_results
    
    end type Hard_Spheres_Observables
    
    type, extends(Hard_Spheres_Observables), public :: Dipolar_Hard_Spheres_Observables
        
        ! Rotate
        integer :: rotate_num_hits = 0
        integer :: rotate_Nreject = 0
        real(DP) :: rotate_reject = 0._DP
        real(DP) :: rotate_rejectSum = 0._DP
        real(DP) :: rotate_rejectAdapt = 0._DP
        real(DP) :: rotate_rejectAvg = 0._DP
        
    end type Dipolar_Hard_Spheres_Observables
    
contains
    
    !> Update rejection
    
    pure subroutine Hard_Spheres_Observables_update_rejections(this)
    
        class(Hard_Spheres_Observables), intent(inout) :: this
        
        this%move_reject = real(this%move_Nreject, DP)/real(this%move_num_hits, DP)
        this%move_Nreject = 0
        this%move_num_hits = 0
        
         select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                this%rotate_reject = real(this%rotate_Nreject, DP)/real(this%rotate_num_hits, DP)
                this%rotate_Nreject = 0
                this%rotate_num_hits = 0
        end select
    
    end subroutine Hard_Spheres_Observables_update_rejections
    
    subroutine Hard_Spheres_Observables_write(this, i_step, obs_unit)
    
        class(Hard_Spheres_Observables), intent(in) :: this
        integer, intent(in) :: i_step, obs_unit
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(obs_unit, *) i_step, this%potential, this%activ, this%move_reject, &
                                   this%rotate_reject
            class default
                write(obs_unit, *) i_step, this%potential, this%activ, this%move_reject
        end select
    
    end subroutine Hard_Spheres_Observables_write
    
    !> Accumulations
    
    pure subroutine Hard_Spheres_Observables_accumulate(this)
    
        class(Hard_Spheres_Observables), intent(inout) :: this
        
        this%potential_sum = this%potential_sum + this%potential
        this%activSum = this%activSum + this%activ
        this%move_rejectSum = this%move_rejectSum + this%move_reject
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                this%rotate_rejectSum = this%rotate_rejectSum + this%rotate_reject
        end select
    
    end subroutine Hard_Spheres_Observables_accumulate
    
    !> Results
    
    subroutine Hard_Spheres_Observables_write_results(this, temperature, num_equilibrium_steps, &
                                                      report_unit)

        class(Hard_Spheres_Observables), intent(in) :: this
        real(DP), intent(in) :: temperature
        integer, intent(in) :: num_equilibrium_steps
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiEx
            
        write(report_unit, *) "Results: "
        
        write(report_unit, *) "    average energy = ", this%potential_sum/real(num_equilibrium_steps, DP)
        potChiEx = -temperature*log(this%activSum/real(num_equilibrium_steps, DP))
        write(report_unit, *) "    average excess chemical potential = ", potChiEx
        write(report_unit, *) "    move rejection rate = ", &
                                   this%move_rejectSum/real(num_equilibrium_steps, DP)
        
        select type (this)
            type is (Dipolar_Hard_Spheres_Observables)
                write(report_unit, *) "    rotate rejection rate = ", &
                                      this%rotate_rejectSum/real(num_equilibrium_steps, DP)
        end select
    
    end subroutine Hard_Spheres_Observables_write_results

end module class_hard_spheres_observables
