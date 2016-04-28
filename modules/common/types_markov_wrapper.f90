module types_markov_wrapper

use types_changes_wrapper, only: Changes_Wrapper
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper

implicit none

private

    type, public :: Markov_Wrapper
        type(Changes_Wrapper) :: changes
        type(Metropolis_Algorithms_Wrapper) :: metropolis_algorithms
    end type Markov_Wrapper

end module types_markov_wrapper
