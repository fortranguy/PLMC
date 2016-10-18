module types_exploring_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Change_Counter
use types_observables_energies, only: Concrete_Observables_Energies

implicit none

private

    !> @note beta_pressures_excess becomes the lateral excess pressures in 2D
    type, public :: Exploring_Observables_Wrapper
        real(DP), allocatable :: maximum_boxes_compression_delta(:)
        real(DP), allocatable :: beta_pressures_excess(:) !> \( \beta p_\text{ex} \)
        real(DP), allocatable :: inv_pow_activities(:, :) !> \( a^{-N} \)
        type(Concrete_Observables_Energies), allocatable :: energies(:)
        type(Concrete_Change_Counter), allocatable :: insertion_counters(:, :)
        real(DP), allocatable :: insertion_successes(:, :)
    end type Exploring_Observables_Wrapper

end module types_exploring_observables_wrapper
