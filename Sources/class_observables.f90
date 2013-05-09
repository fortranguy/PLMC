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
        real(DP) :: rej
        real(DP) :: rejSum
        real(DP) :: rejAdapt
        
        ! Potential energy
        real(DP) :: Epot
        real(DP) :: EpotSum
        
        ! Inverse of activity
        real(DP) :: activ
        real(DP) :: activSum
    
    contains
    
        procedure :: init => Observables_init
        procedure :: results => Observables_results
    
    end type Observables
    
    type, extends(Observables), public :: MoreObservables
        
        ! Rotate
        integer :: Nrotate
        
        ! Rejection
        integer :: NrejRot
        real(DP) :: rejRot
        real(DP) :: rejRotSum
        real(DP) :: rejRotAdapt
        
    end type MoreObservables
    
contains

    subroutine Observables_init(this)
        
        class(Observables), intent(out) :: this
        
        this%Nmove = 0
        
        this%Nrej = 0
        this%rej = 0._DP
        this%rejSum = 0._DP
        this%rejAdapt = 0._DP
        
        this%EpotSum = 0._DP        
        this%activSum = 0._DP
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
                
                this%Nrotate = 0
                
                this%NrejRot = 0
                this%rejRot = 0._DP
                this%rejRotSum = 0._DP
                this%rejRotAdapt = 0._DP
                
        end select
        
    end subroutine Observables_init
    
    !> Results
    
    subroutine Observables_results(this, Ncol, report_unit)

        class(Observables), intent(in) :: this
        integer, intent(in) :: Ncol
        integer, intent(in) :: report_unit
        
        real(DP) :: potChiId, potChiEx
            
        write(report_unit, *) "Results :"
        
        write(report_unit, *) "    average energy = ", this%EpotSum/real(Nstep, DP)
        write(report_unit, *) "    average energy per particule = ", &
                                   this%EpotSum/real(Nstep, DP)/real(Ncol, DP)
            
        potChiId = -Tstar*log( product(Lsize)/real(Ncol+1,DP) )
        write(report_unit, *) "    ideal chemical potential = ", potChiId
        potChiEx = -Tstar*log( this%activSum/real(Nstep, DP) )
        write(report_unit, *) "    average excess chemical potential = ", potChiEx           
        write(report_unit, *) "    potChi.avg = ", potChiId + potChiEx
        
        write(report_unit, *) "    Rejection rate = ", this%rejSum/real(Nstep, DP)
        
        select type (this)
        
            type is (Observables)
            
            class is (MoreObservables)
            
                write(report_unit, *) "not yet implemented"
            
        end select
    
    end subroutine Observables_results

end module class_observables
