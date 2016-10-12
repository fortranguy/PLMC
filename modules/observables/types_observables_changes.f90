module types_observables_changes

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line
use module_changes_success, only: Concrete_Change_Counter, Concrete_Changes_Counter, &
    Concrete_Changes_Success, Concrete_Change_Counters_Line

implicit none

private

    type, public :: Concrete_Observables_Changes
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Change_Counters_Line), allocatable :: switches_counters(:)
        type(Reals_Line), allocatable :: switches_successes(:)
        type(Concrete_Change_Counter), allocatable :: transmutations_counters(:, :)
        real(DP), allocatable :: transmutations_successes(:, :)
    end type Concrete_Observables_Changes

end module types_observables_changes
