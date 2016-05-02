module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success, &
    Concrete_Switch_Counters
use types_reals_line, only: Reals_Line

implicit none

private

    type, public :: Generating_Observables_Wrapper
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Switch_Counters), allocatable :: switches_counters(:)
        type(Reals_Line), allocatable :: switches_successes(:)
        real(DP), allocatable :: field_energies(:), walls_energies(:)
        type(Reals_Line), allocatable :: short_energies(:), dipolar_energies(:)
        real(DP) :: dipolar_mixture_energy = 0._DP
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
