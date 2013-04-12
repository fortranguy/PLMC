!> \brief Description of the Observables class

module class_observables

use data_constants
use data_mc

implicit none

private

    type, public :: Observables
    
        !Rejection
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
        procedure :: add => Observables_add
        procedure :: addReject => Observables_addReject
        procedure :: results => Observables_results
    
    end type Observables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nrejects = 0
        this%rejectsRateSum = 0._DP
        
        this%ePot_totalSum = 0._DP        
        this%activExInvSum = 0._DP
        
    end subroutine Observables_init
    
    subroutine Observables_add(this)
    
        class(Observables), intent(inout) :: this
    
        this%ePot_totalSum = this%ePot_totalSum + this%ePot_total
        this%activExInvSum = this%activExInvSum + this%activExInv
            
    end subroutine Observables_add
    
    subroutine Observables_addReject(this)
    
        class(Observables), intent(inout) :: this
        
        this%rejectsRateSum = this%rejectsRateSum + real(this%Nrejects, DP) / &
            real(Nmove, DP)
        this%Nrejects = 0
    
    end subroutine Observables_addReject
    
    !> Results
    
    subroutine Observables_results(this, duration, unitRapport)

        class(Observables), intent(inout) :: this
        real(DP), intent(in) :: duration
        integer, intent(in) :: unitRapport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitRapport, *) "Results :"
        
        write(unitRapport, *) "    average energy = ", &
            this%ePot_totalSum/realNstep
        write(unitRapport, *) "    average energy per particule = ", &
            this%ePot_totalSum/realNstep/real(inter_Ncol, DP)
            
        potChiId = -Tstar*log( product(Lsize)/real(inter_Ncol+1,DP) )
        write(unitRapport, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( this%activExInvSum/realNstep )
        write(unitRapport, *) "    average excess chemical potential = ", &
            potChiEx           
        write(unitRapport, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(unitRapport, *) "    Rejection rate = ", &
            this%rejectsRateSum/real(Nstep+Ntherm, DP)
        write(unitRapport, *) "    duration =", duration/60._DP, "min"        
    
    end subroutine Observables_results

end module class_observables
