!> \brief Description of the Observables class

module class_observables

use data_precisions, only: DP

implicit none
private

    type, public :: Observables
    
        ! Move
        integer :: move_Nhit = 0
        integer :: move_Nreject = 0
        real(DP) :: move_reject = 0._DP
        real(DP) :: move_rejectSum = 0._DP
        real(DP) :: move_rejectAdapt = 0._DP
        real(DP) :: move_rejectAvg = 0._DP
        
        ! Potential energy
        real(DP) :: Epot
        real(DP) :: EpotSum = 0._DP
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum = 0._DP
    
    contains
    
        procedure :: update_rejections => Observables_update_rejections
        procedure :: write => Observables_write
        procedure :: accumulate => Observables_accumulate
        procedure :: write_results => Observables_write_results
    
    end type Observables
    
    type, extends(Observables), public :: MoreObservables
        
        ! Rotate
        integer :: rotate_Nhit = 0
        integer :: rotate_Nreject = 0
        real(DP) :: rotate_reject = 0._DP
        real(DP) :: rotate_rejectSum = 0._DP
        real(DP) :: rotate_rejectAdapt = 0._DP
        real(DP) :: rotate_rejectAvg = 0._DP
        
    end type MoreObservables
    
contains
    
    !> Update rejection
    
    pure subroutine Observables_update_rejections(this)
    
        class(Observables), intent(inout) :: this
        
        this%move_reject = real(this%move_Nreject, DP)/real(this%move_Nhit, DP)
        this%move_Nreject = 0
        this%move_Nhit = 0
        
         select type (this)
            type is (MoreObservables)
                this%rotate_reject = real(this%rotate_Nreject, DP)/real(this%rotate_Nhit, DP)
                this%rotate_Nreject = 0
                this%rotate_Nhit = 0
        end select
    
    end subroutine Observables_update_rejections
    
    subroutine Observables_write(this, iStep, obs_unit)
    
        class(Observables), intent(in) :: this
        integer, intent(in) :: iStep, obs_unit
        
        select type (this)
            type is (MoreObservables)
                write(obs_unit, *) iStep, this%Epot, this%activ, this%move_reject, &
                                   this%rotate_reject
            class default
                write(obs_unit, *) iStep, this%Epot, this%activ, this%move_reject
        end select
    
    end subroutine Observables_write
    
    !> Accumulations
    
    pure subroutine Observables_accumulate(this)
    
        class(Observables), intent(inout) :: this
        
        this%EpotSum = this%EpotSum + this%Epot
        this%activSum = this%activSum + this%activ
        this%move_rejectSum = this%move_rejectSum + this%move_reject
        
        select type (this)
            type is (MoreObservables)
                this%rotate_rejectSum = this%rotate_rejectSum + this%rotate_reject
        end select
    
    end subroutine Observables_accumulate
    
    !> Results
    
    subroutine Observables_write_results(this, temperature, num_equilibrium_steps, report_unit)

        class(Observables), intent(in) :: this
        real(DP), intent(in) :: temperature
        integer, intent(in) :: num_equilibrium_steps
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiEx
            
        write(report_unit, *) "Results: "
        
        write(report_unit, *) "    average energy = ", this%EpotSum/real(num_equilibrium_steps, DP)
        potChiEx = -temperature*log(this%activSum/real(num_equilibrium_steps, DP))
        write(report_unit, *) "    average excess chemical potential = ", potChiEx
        write(report_unit, *) "    move rejection rate = ", &
                                   this%move_rejectSum/real(num_equilibrium_steps, DP)
        
        select type (this)
            type is (MoreObservables)
                write(report_unit, *) "    rotate rejection rate = ", &
                                      this%rotate_rejectSum/real(num_equilibrium_steps, DP)
        end select
    
    end subroutine Observables_write_results

end module class_observables
