module types_exploring_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
implicit none

private

    type, public :: Exploring_Observables_Wrapper
        real(DP), allocatable :: inv_pow_activities(:)
    end type Exploring_Observables_Wrapper

end module types_exploring_observables_wrapper
