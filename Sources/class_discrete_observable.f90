module class_discrete_observable

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    ! Observables

    type, public :: Discrete_Observables
    
        integer :: num_hits = 0
        integer :: num_rejections = 0
        real(DP) :: rejection_rate = 0._DP
        real(DP) :: sum_rejection = 0._DP
        
    contains
    
        procedure :: update_rejection => Discrete_Observables_update_rejection
        
    end type Discrete_Observables
    
    type, extends(Discrete_Observables), public :: Adapting_Discrete_Observables
    
        real(DP) :: rejection_adapt = 0._DP
        real(DP) :: rejection_average = 0._DP
        
    contains
    
        procedure :: accumulate_rejection => Adapting_Discrete_Observables_accumulate_rejection
        procedure :: average_rejection => Adapting_Discrete_Observables_average_rejection
        
    end type Adapting_Discrete_Observables
    
contains

    subroutine Discrete_Observables_update_rejection(this)
    
        class(Discrete_Observables), intent(inout) :: this    

        this%rejection_rate = real(this%num_rejections, DP) / real(this%num_hits, DP)
        this%num_rejections = 0
        this%num_hits = 0
        
    end subroutine Discrete_Observables_update_rejection
    
    subroutine Adapting_Discrete_Observables_accumulate_rejection(this)
    
        class(Adapting_Discrete_Observables), intent(inout) :: this   
    
        this%rejection_adapt = this%rejection_adapt + this%rejection_rate
        
    end subroutine Adapting_Discrete_Observables_accumulate_rejection
    
    subroutine Adapting_Discrete_Observables_average_rejection(this, period_adaptation)
    
        class(Adapting_Discrete_Observables), intent(inout) :: this
        integer, intent(in) :: period_adaptation
    
        this%rejection_average = this%rejection_adapt / real(period_adaptation - 1, DP)
        this%rejection_adapt = 0._DP
        
    end subroutine Adapting_Discrete_Observables_average_rejection
    
end module class_discrete_observable
