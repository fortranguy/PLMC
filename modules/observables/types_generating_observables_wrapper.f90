module types_generating_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_reals_line, only: Reals_Line
use module_changes_success, only: Concrete_Change_Counter, Concrete_Changes_Counter, &
    Concrete_Changes_Success, Concrete_Change_Counters_Line
use types_observables_energies, only: Concrete_Observables_Energies
use types_observables_changes, only: Concrete_Observables_Changes

implicit none

private

    !> @todo
    !> second part to remove
    type, public :: Generating_Observables_Wrapper
        real(DP), allocatable :: accessible_domains_size(:, :)
        type(Concrete_Change_Counters_Line), allocatable :: volumes_change_counter(:)
        type(Reals_Line), allocatable :: volumes_change_success(:)
        integer, allocatable :: gemc_nums_particles(:, :)
        type(Concrete_Observables_Energies), allocatable :: gemc_energies(:)
        type(Concrete_Observables_Changes), allocatable :: changes(:)

        real(DP) :: accessible_domain_size(num_dimensions) = 0._DP
        type(Concrete_Change_Counter) :: volume_change_counter
        real(DP) :: volume_change_success = 0._DP
        integer, allocatable :: nums_particles(:)
        type(Concrete_Observables_Energies) :: energies
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Change_Counters_Line), allocatable :: switches_counters(:)
        type(Reals_Line), allocatable :: switches_successes(:)
        type(Concrete_Change_Counter), allocatable :: transmutations_counters(:, :)
        real(DP), allocatable :: transmutations_successes(:, :)
    end type Generating_Observables_Wrapper

end module types_generating_observables_wrapper
