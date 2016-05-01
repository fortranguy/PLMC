module types_markov_chain_generator_wrapper

use types_changes_wrapper, only: Changes_Wrapper
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator

implicit none

private

    type, public :: Markov_Chain_Generator_Wrapper
        type(Changes_Wrapper) :: changes
        type(Metropolis_Algorithms_Wrapper) :: metropolis_algorithms
        class(Abstract_PLMC_Propagator), allocatable :: plmc_propagator
    end type Markov_Chain_Generator_Wrapper

end module types_markov_chain_generator_wrapper
