module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use module_changes_success, only: Concrete_Change_Counter, Concrete_Changes_Counter, &
    Concrete_Changes_Success, Concrete_Change_Counters_Line
use types_reals_line, only: Reals_Line
use types_observables_energies, only: Concrete_Energies

implicit none

private

    type, public :: Generating_Observables_Wrapper
        real(DP) :: accessible_domain_size(num_dimensions) = 0._DP
        type(Concrete_Change_Counter) :: volume_change_counter
        real(DP) :: volume_change_success = 0._DP
        integer, allocatable :: nums_particles(:)
        type(Concrete_Energies) :: energies
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Change_Counters_Line), allocatable :: switches_counters(:)
        type(Reals_Line), allocatable :: switches_successes(:)
        type(Concrete_Change_Counter), allocatable :: transmutations_counters(:, :)
        real(DP), allocatable :: transmutations_successes(:, :)
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
