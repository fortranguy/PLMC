module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_real_line, only: Real_Line
use module_changes_success, only: Concrete_Change_Counter, Concrete_Change_Counter_Line
use types_observables_energies, only: Concrete_Observables_Energies
use types_observables_changes, only: Concrete_Observables_Changes

implicit none

private

    type, public :: Generating_Observables_Wrapper
        real(DP), allocatable :: accessible_domains_size(:, :)
        type(Concrete_Change_Counter), allocatable :: teleportations_counters(:, :, :)
        real(DP), allocatable :: teleportations_successes(:, :, :)
        type(Concrete_Change_Counter_Line), allocatable :: volumes_change_counter(:)
        type(Real_Line), allocatable :: volumes_change_success(:)
        integer, allocatable :: nums_particles(:, :)
        type(Concrete_Observables_Energies), allocatable :: energies(:)
        type(Concrete_Observables_Changes), allocatable :: changes(:)
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
