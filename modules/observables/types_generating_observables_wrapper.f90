module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line
use module_changes_success, only: Concrete_Change_Counters_Line
use types_observables_energies, only: Concrete_Observables_Energies
use types_observables_changes, only: Concrete_Observables_Changes

implicit none

private

    type, public :: Generating_Observables_Wrapper
        real(DP), allocatable :: accessible_domains_size(:, :)
        type(Concrete_Change_Counters_Line), allocatable :: volumes_change_counter(:)
        type(Reals_Line), allocatable :: volumes_change_success(:)
        integer, allocatable :: nums_particles(:, :)
        type(Concrete_Observables_Energies), allocatable :: energies(:)
        type(Concrete_Observables_Changes), allocatable :: changes(:)
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
