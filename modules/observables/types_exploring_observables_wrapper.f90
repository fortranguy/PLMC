module types_exploring_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Change_Counter
use types_reals_line, only: Reals_Line
use types_observables_energies, only: Concrete_Energies

implicit none

private

    !> beta_pressure_excess becomes the lateral excess pressure in 2D
    type, public :: Exploring_Observables_Wrapper
        real(DP) :: maximum_box_compression_delta = 0._DP
        real(DP) :: beta_pressure_excess = 0._DP !! \( \beta p_\text{ex} \)
        real(DP), allocatable :: inv_pow_activities(:) !! \( a^{-N} \)
        type(Concrete_Energies) :: energies
        type(Concrete_Change_Counter), allocatable :: insertion_counters(:)
        real(DP), allocatable :: insertion_successes(:)
    end type Exploring_Observables_Wrapper

end module types_exploring_observables_wrapper
