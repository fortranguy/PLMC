module classes_plmc_propagator

use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithm_Pointer, &
    Metropolis_Algorithms_Wrapper

implicit none

private

    type, abstract, public :: Abstract_PLMC_Propagator
    private
    contains
    end type Abstract_PLMC_Propagator

end module classes_plmc_propagator
