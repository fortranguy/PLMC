module types_mixture_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use module_changes_success, only:Concrete_Mixture_Changes_Counters, Concrete_Mixture_Changes_Success
use module_particles_energy, only: Concrete_Particles_Energy

implicit none

private

    type, public :: Concrete_Mixture_Observables
        type(Concrete_Mixture_Changes_Counters) :: changes_counters
        type(Concrete_Mixture_Changes_Success) :: changes_success
        type(Concrete_Particles_Energy) :: particles_energies(2)
        real(DP) :: inter_energy
    end type Concrete_Mixture_Observables

end module types_mixture_observables
