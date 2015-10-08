module class_candidates_selector

use data_constants, only: num_components
use procedures_random, only: random_integer
implicit none

private

    type, abstract, public :: Abstract_Candidates_Selector
    contains
        procedure, nopass :: set => Abstract_Candidates_Selector_set
    end type Abstract_Candidates_Selector

contains

    subroutine Abstract_Candidates_Selector_set(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = random_integer(num_components)
        i_spectator = mod(i_actor, num_components) + 1
    end subroutine Abstract_Candidates_Selector_set

end module class_candidates_selector

