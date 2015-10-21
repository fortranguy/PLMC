module types_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success
use module_component_energy, only: Concrete_Component_Energy, Concrete_Inter_Energy

implicit none

private

    type, public :: Observables_Wrapper_old
        type(Concrete_Changes_Counter) :: changes_counter
        type(Concrete_Changes_Success) :: changes_success
        type(Concrete_Component_Energy) :: component_energy ! to remove
    end type Observables_Wrapper_old

    type, public :: Concrete_Inter_Energies
        real(DP), allocatable :: with_components(:)
    end type Concrete_Inter_Energies

    type, public :: Observables_Wrapper
        type(Observables_Wrapper_old) :: intras(num_components) ! to remove
        type(Concrete_Inter_Energy) :: inter_energy ! to remove
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Inter_Energies), allocatable :: short_energies(:)
        real(DP), allocatable :: walls_energies(:)
        type(Concrete_Inter_Energies), allocatable :: long_energies(:)
        real(DP), allocatable :: field_energies(:)
    end type Observables_Wrapper

end module types_observables_wrapper
