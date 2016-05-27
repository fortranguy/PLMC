module types_exploring_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Change_Counter

implicit none

private

    type, public :: Exploring_Observables_Wrapper
        type(Concrete_Change_Counter), allocatable :: widom_counters(:)
        real(DP), allocatable :: widom_successes(:)
        real(DP), allocatable :: inv_pow_activities(:)
    end type Exploring_Observables_Wrapper

end module types_exploring_observables_wrapper
