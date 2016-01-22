module types_observables_wrapper

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success

implicit none

private

    type, public :: Concrete_Components_Energies
        real(DP), allocatable :: with_components(:)
    end type Concrete_Components_Energies

    type, public :: Observables_Wrapper
        type(Concrete_Changes_Counter), allocatable :: changes_counters(:)
        type(Concrete_Changes_Success), allocatable :: changes_sucesses(:)
        type(Concrete_Components_Energies), allocatable :: short_energies(:)
        real(DP), allocatable :: walls_energies(:)
        type(Concrete_Components_Energies), allocatable :: long_energies(:)
        real(DP) :: long_mixture_energy
        real(DP), allocatable :: field_energies(:)
    end type Observables_Wrapper

end module types_observables_wrapper
