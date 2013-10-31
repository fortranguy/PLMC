!> \brief Description of the Observables class

module class_observables

use data_precisions, only : DP
use data_box, only : Volume
use data_monteCarlo, only : Temperature, Nstep

implicit none
private

    type, public :: Observables
    
        ! Move
        integer :: Nmove
        integer :: move_Nreject
        real(DP) :: move_reject
        real(DP) :: move_rejectSum
        real(DP) :: move_rejectAdapt
        real(DP) :: move_rejectAvg
        
        ! Potential energy
        real(DP) :: Epot
        real(DP) :: EpotSum
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum
    
    contains
    
        procedure :: init => Observables_init
        procedure :: update_rejections => Observables_update_rejections
        procedure :: accumulate => Observables_accumulate
        procedure :: print_results => Observables_print_results
    
    end type Observables
    
    type, extends(Observables), public :: MoreObservables
        
        ! Rotate
        integer :: Nrotate
        integer :: rotate_Nreject
        real(DP) :: rotate_reject
        real(DP) :: rotate_rejectSum
        real(DP) :: rotate_rejectAdapt
        real(DP) :: rotate_rejectAvg
        
    end type MoreObservables
    
contains

    !> Initialisation

    pure subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nmove = 0        
        this%move_Nreject = 0
        this%move_reject = 0._DP
        this%move_rejectSum = 0._DP
        this%move_rejectAdapt = 0._DP
        this%move_rejectAvg = 0._DP
        
        this%EpotSum = 0._DP        
        this%activSum = 0._DP
        
        select type (this)
            
            type is (MoreObservables)
                
                this%Nrotate = 0                
                this%rotate_Nreject = 0
                this%rotate_reject = 0._DP
                this%rotate_rejectSum = 0._DP
                this%rotate_rejectAdapt = 0._DP
                this%rotate_rejectAvg = 0._DP
                
        end select
        
    end subroutine Observables_init
    
    !> Update rejection
    
    pure subroutine Observables_update_rejections(this)
    
        class(Observables), intent(inout) :: this
        
        this%move_reject = real(this%move_Nreject, DP)/real(this%Nmove, DP)
        this%move_Nreject = 0
        this%Nmove = 0
        
         select type (this)
        
            type is (MoreObservables)
                
                this%rotate_reject = real(this%rotate_Nreject, DP)/real(this%Nrotate, DP)
                this%rotate_Nreject = 0
                this%Nrotate = 0
        
        end select
    
    end subroutine Observables_update_rejections
    
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
    
    subroutine Observables_print_results(this, report_unit)

        class(Observables), intent(in) :: this
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiEx
            
        write(report_unit, *) "Results :"
        
        write(report_unit, *) "    average energy = ", this%EpotSum/real(Nstep, DP)
            
        potChiEx = -Temperature*log(this%activSum/real(Nstep, DP))
        write(report_unit, *) "    average excess chemical potential = ", potChiEx
        
        write(report_unit, *) "    move rejection rate = ", this%move_rejectSum/real(Nstep, DP)
        
        select type (this)
            
            type is (MoreObservables)
            
                write(report_unit, *) "    rotate rejection rate = ", &
                                      this%rotate_rejectSum/real(Nstep, DP)
            
        end select
    
    end subroutine Observables_print_results

end module class_observables
