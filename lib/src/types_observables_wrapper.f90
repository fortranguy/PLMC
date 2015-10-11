module types_observables_wrapper

use data_constants, only: num_components
use module_changes_success, only: Concrete_Changes_Counter, Concrete_Changes_Success
use module_particles_energy, only: Concrete_Particles_Energy, Concrete_Inter_Energy

implicit none

private

    type, public :: Observables_Wrapper
        type(Concrete_Changes_Counter) :: changes_counter
        type(Concrete_Changes_Success) :: changes_success
        type(Concrete_Particles_Energy) :: particles_energy
    end type Observables_Wrapper

    type, public :: Mixture_Observables_Wrapper
        type(Observables_Wrapper) :: intras(num_components)
        type(Concrete_Inter_Energy) :: inter_energy
    end type Mixture_Observables_Wrapper

end module types_observables_wrapper
