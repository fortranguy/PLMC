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
    
        ! Move Rejection
        integer :: Nmove_reject
        real(DP) :: move_reject
        real(DP) :: move_rejectSum
        real(DP) :: move_rejectAdapt
        
        ! Potential energy
        real(DP) :: Epot
        real(DP) :: EpotSum
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum
    
    contains
    
        procedure :: init => Observables_init
        procedure :: printResults => Observables_printResults
    
    end type Observables
    
    type, extends(Observables), public :: MoreObservables
        
        ! Rotate
        integer :: Nrotate
        
        ! Rotate rejection
        integer :: Nrotate_reject
        real(DP) :: rotate_reject
        real(DP) :: rotate_rejectSum
        real(DP) :: rotate_rejectAdapt
        
    end type MoreObservables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nmove = 0
        
        this%Nmove_reject = 0
        this%move_reject = 0._DP
        this%move_rejectSum = 0._DP
        this%move_rejectAdapt = 0._DP
        
        this%EpotSum = 0._DP        
        this%activSum = 0._DP
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
                
                this%Nrotate = 0
                
                this%Nrotate_reject = 0
                this%rotate_reject = 0._DP
                this%rotate_rejectSum = 0._DP
                this%rotate_rejectAdapt = 0._DP
                
        end select
        
    end subroutine Observables_init
    
    !> Results
    
    subroutine Observables_printResults(this, Ncol, report_unit)

        class(Observables), intent(in) :: this
        integer, intent(in) :: Ncol
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiId, potChiEx
            
        write(report_unit, *) "Results :"
        
        write(report_unit, *) "    average energy = ", this%EpotSum/real(Nstep, DP)
            
        potChiId = -Temperature*log(Volume/real(Ncol+1,DP))
        write(report_unit, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Temperature*log(this%activSum/real(Nstep, DP))
        write(report_unit, *) "    average excess chemical potential = ", potChiEx           
        write(report_unit, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(report_unit, *) "    move_rejection rate = ", this%move_rejectSum/real(Nstep, DP)
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
            
                write(report_unit, *) "    Rotation move_rejection rate = ", &
                                      this%rotate_rejectSum/real(Nstep, DP)
            
        end select
    
    end subroutine Observables_printResults

end module class_observables
