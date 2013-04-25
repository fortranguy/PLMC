!> \brief Description of the Observables class

module class_observables

use data_constants
use data_mc

implicit none

private

    type, public :: Observables
    
        ! Move
        integer :: Nmove
    
        ! Rejection
        integer :: Nrej
        real(DP) :: rejRateSum
        
        ! Potential energy
        real(DP) :: ePot
        real(DP) :: ePotSum
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum
    
    contains
    
        procedure :: init => Observables_init
        procedure :: addPhysical => Observables_addPhysical
        procedure :: addReject => Observables_addReject
        procedure :: results => Observables_results
    
    end type Observables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nmove = 0
        
        this%Nrej = 0
        this%rejRateSum = 0._DP
        
        this%ePotSum = 0._DP        
        this%activSum = 0._DP
        
    end subroutine Observables_init
    
    subroutine Observables_addPhysical(this)
    
        class(Observables), intent(inout) :: this
    
        this%ePotSum = this%ePotSum + this%ePot
        this%activSum = this%activSum + this%activ
            
    end subroutine Observables_addPhysical
    
    subroutine Observables_addReject(this)
    
        class(Observables), intent(inout) :: this
        
        this%rejRateSum = this%rejRateSum + real(this%Nrej, DP) / &
            real(this%Nmove, DP)
        this%Nrej = 0
        this%Nmove = 0
    
    end subroutine Observables_addReject
    
    !> Results
    
    subroutine Observables_results(this, Ncol, unitReport)

        class(Observables), intent(in) :: this
        integer, intent(in) :: Ncol
        integer, intent(in) :: unitReport
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(unitReport, *) "Results :"
        
        write(unitReport, *) "    average energy = ", &
            this%ePotSum/realNstep
        write(unitReport, *) "    average energy per particule = ", &
            this%ePotSum/realNstep/real(Ncol, DP)
            
        potChiId = -Tstar*log( product(Lsize)/real(Ncol+1,DP) )
        write(unitReport, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( this%activSum/realNstep )
        write(unitReport, *) "    average excess chemical potential = ", &
            potChiEx           
        write(unitReport, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(unitReport, *) "    Rejection rate = ", &
            this%rejRateSum/real(Nstep+Ntherm, DP)
    
    end subroutine Observables_results

end module class_observables
