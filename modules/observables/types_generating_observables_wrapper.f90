module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_real_wrapper, only: Real_Line, Real_Triangle_Line
use module_changes_success, only: Concrete_Change_Counter, Concrete_Change_Counter_Line, &
    Concrete_Change_Counter_Triangle_Line
use types_observables_energies, only: Concrete_Observables_Energies
use types_observables_changes, only: Concrete_Observables_Changes

implicit none

private

    type, public :: Generating_Observables_Wrapper
        real(DP), allocatable :: accessible_domains_size(:, :)
        type(Concrete_Change_Counter), allocatable :: volumes_change_counter(:)
        real(DP), allocatable :: volumes_change_success(:)
        type(Concrete_Change_Counter_Line), allocatable :: volumes_exchange_counter(:)
        type(Real_Line), allocatable :: volumes_exchange_success(:)
        type(Concrete_Change_Counter), allocatable :: teleportations_counters(:, :, :)
        real(DP), allocatable :: teleportations_successes(:, :, :)
        type(Concrete_Change_Counter_Triangle_Line), allocatable :: switches_counters(:)
        type(Real_Triangle_Line), allocatable :: switches_successes(:)
        integer, allocatable :: nums_particles(:, :)
        type(Concrete_Observables_Energies), allocatable :: energies(:)
        type(Concrete_Observables_Changes), allocatable :: changes(:)
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
