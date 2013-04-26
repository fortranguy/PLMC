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
    
    !> Results
    
    subroutine Observables_results(this, Ncol, report_unit)

        class(Observables), intent(in) :: this
        integer, intent(in) :: Ncol
        integer, intent(in) :: report_unit
        
        real(DP) :: realNstep = real(Nstep, DP)
        real(DP) :: potChiId, potChiEx
    
        write(report_unit, *) "Results :"
        
        write(report_unit, *) "    average energy = ", &
            this%ePotSum/realNstep
        write(report_unit, *) "    average energy per particule = ", &
            this%ePotSum/realNstep/real(Ncol, DP)
            
        potChiId = -Tstar*log( product(Lsize)/real(Ncol+1,DP) )
        write(report_unit, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( this%activSum/realNstep )
        write(report_unit, *) "    average excess chemical potential = ", &
            potChiEx           
        write(report_unit, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(report_unit, *) "    Rejection rate = ", &
            this%rejRateSum/real(Nstep+Ntherm, DP)
    
    end subroutine Observables_results

end module class_observables
