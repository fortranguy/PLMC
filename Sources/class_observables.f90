!> \brief Description of the Observables class

module class_observables

use data_precisions, only : DP
use data_constants
use data_cell, only : Volume
use data_mc

implicit none

private

    type, public :: Observables
    
        ! Move
        integer :: Nmove
    
        ! Rejection
        integer :: Nreject
        real(DP) :: reject
        real(DP) :: rejectSum
        real(DP) :: rejectAdapt
        
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
        
        ! Rejection
        integer :: NrejectRot
        real(DP) :: rejectRot
        real(DP) :: rejectRotSum
        real(DP) :: rejectRotAdapt
        
    end type MoreObservables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nmove = 0
        
        this%Nreject = 0
        this%reject = 0._DP
        this%rejectSum = 0._DP
        this%rejectAdapt = 0._DP
        
        this%EpotSum = 0._DP        
        this%activSum = 0._DP
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
                
                this%Nrotate = 0
                
                this%NrejectRot = 0
                this%rejectRot = 0._DP
                this%rejectRotSum = 0._DP
                this%rejectRotAdapt = 0._DP
                
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
        
        write(report_unit, *) "    Rejection rate = ", this%rejectSum/real(Nstep, DP)
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
            
                write(report_unit, *) "    Rotation rejection rate = ", &
                                      this%rejectRotSum/real(Nstep, DP)
            
        end select
    
    end subroutine Observables_printResults

end module class_observables
