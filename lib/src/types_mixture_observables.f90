module types_mixture_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_change_counter, only: Concrete_Change_Counter
use module_particles_energy, only: Concrete_Particles_Energy

implicit none

private

    type, public :: Concrete_Mixture_Observables
        type(Concrete_Change_Counter) :: move_counters(2)
        real(DP) :: move_success_ratio(2)
        type(Concrete_Particles_Energy) :: particles_energies(2)
        real(DP) :: inter_energy
    end type Concrete_Mixture_Observables

end module types_mixture_observables
