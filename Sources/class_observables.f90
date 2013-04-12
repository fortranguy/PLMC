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
        procedure :: addRej => Observables_addRej
    
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
    
    subroutine Observables_addRej(this)
    
        class(Observables), intent(inout) :: this
        
        this%rejectsRateSum = this%rejectsRateSum + real(this%Nrejects, DP) / &
            real(Nmove, DP)
        this%Nrejects = 0
    
    end subroutine Observables_addRej

end module class_observables
