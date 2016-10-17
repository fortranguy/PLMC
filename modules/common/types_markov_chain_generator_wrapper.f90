module types_markov_chain_generator_wrapper

use types_changes_wrapper, only: Changes_Wrapper
use classes_generating_algorithm, only: Generating_Algorithm_Wrapper
use types_generating_algorithms_wrapper, only: Generating_Algorithms_Wrapper
use classes_plmc_propagator, only: Abstract_PLMC_Propagator

implicit none

private

    type, public :: Markov_Chain_Generator_Wrapper
        type(Changes_Wrapper) :: changes
        type(Generating_Algorithm_Wrapper), allocatable :: gemc_generating_algorithms(:)
        class(Abstract_PLMC_Propagator), allocatable :: plmc_propagator

        type(Generating_Algorithms_Wrapper) :: generating_algorithms
    end type Markov_Chain_Generator_Wrapper

end module types_markov_chain_generator_wrapper
