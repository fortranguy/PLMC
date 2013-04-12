!> \brief Description of the Observables class

module class_observables

use data_constants

implicit none

private

    type, public :: Obversables
    
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
    
    end type Obversables

end module class_observables
