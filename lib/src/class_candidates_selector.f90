module class_candidates_selector

use data_constants, only: num_components
use procedures_random, only: random_integer
implicit none

private

    type, abstract, public :: Abstract_Candidates_Selector
    contains
        procedure, nopass :: set => Abstract_Candidates_Selector_set
    end type Abstract_Candidates_Selector

    type, extends(Abstract_Candidates_Selector), public :: First_Candidate_Selector
    contains
        procedure, nopass :: set => First_Candidate_Selector_set
    end type First_Candidate_Selector

    type, extends(Abstract_Candidates_Selector), public :: Second_Candidate_Selector
    contains
        procedure, nopass :: set => Second_Candidate_Selector_set
    end type Second_Candidate_Selector


    type, extends(Abstract_Candidates_Selector), public :: Null_Candidates_Selector
    contains
        procedure, nopass :: set => Null_Candidates_Selector_set
    end type Null_Candidates_Selector

contains

    subroutine Abstract_Candidates_Selector_set(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = random_integer(num_components)
        i_spectator = mod(i_actor, num_components) + 1
    end subroutine Abstract_Candidates_Selector_set

    subroutine First_Candidate_Selector_set(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 1
        i_spectator = 2
    end subroutine First_Candidate_Selector_set

    subroutine Second_Candidate_Selector_set(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 2
        i_spectator = 1
    end subroutine Second_Candidate_Selector_set

    subroutine Null_Candidates_Selector_set(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator
        i_actor = 0
        i_spectator = 0
    end subroutine Null_Candidates_Selector_set

end module class_candidates_selector

