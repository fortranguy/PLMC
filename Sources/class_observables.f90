!> \brief Description of the Observables class

module class_observables

use data_constants
use data_mc

implicit none

private

    type, public :: Observables
    
        ! Randomly chosen
        
        integer :: NrandMove
    
        ! Rejection
        integer :: Nrejects
        real(DP) :: rejectsRateSum
        
        ! Potential energy
        real(DP) :: ePot_total
        real(DP) :: ePot_totalSum
        
        ! Inverse of activity
        real(DP) :: activExInv
        real(DP) :: activExInvSum
    
    contains
    
        procedure :: init => Observables_init
        procedure :: addPhysical => Observables_addPhysical
        procedure :: addReject => Observables_addReject
        procedure :: results => Observables_results
    
    end type Observables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%NrandMove = 0
        
        this%Nrejects = 0
        this%rejectsRateSum = 0._DP
        
        this%ePot_totalSum = 0._DP        
        this%activExInvSum = 0._DP
        
    end subroutine Observables_init
    
    subroutine Observables_addPhysical(this)
    
        class(Observables), intent(inout) :: this
    
        this%ePot_totalSum = this%ePot_totalSum + this%ePot_total
        this%activExInvSum = this%activExInvSum + this%activExInv
            
    end subroutine Observables_addPhysical
    
    subroutine Observables_addReject(this)
    
        class(Observables), intent(inout) :: this
        
        this%rejectsRateSum = this%rejectsRateSum + real(this%Nrejects, DP) / &
            real(Nmove, DP)
        this%Nrejects = 0
    
    end subroutine Observables_addReject
    
    !> Results
    
    subroutine Observables_results(this, duration, unitReport)

        class(Observables), intent(inout) :: this
        real(DP), intent(in) :: duration
        integer, intent(in) :: unitReport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitReport, *) "Results :"
        
        write(unitReport, *) "    average energy = ", &
            this%ePot_totalSum/realNstep
        write(unitReport, *) "    average energy per particule = ", &
            this%ePot_totalSum/realNstep/real(inter_Ncol, DP)
            
        potChiId = -Tstar*log( product(Lsize)/real(inter_Ncol+1,DP) )
        write(unitReport, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( this%activExInvSum/realNstep )
        write(unitReport, *) "    average excess chemical potential = ", &
            potChiEx           
        write(unitReport, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(unitReport, *) "    Rejection rate = ", &
            this%rejectsRateSum/real(Nstep+Ntherm, DP)
        write(unitReport, *) "    duration =", duration/60._DP, "min"        
    
    end subroutine Observables_results

end module class_observables
