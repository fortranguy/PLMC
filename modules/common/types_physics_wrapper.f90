module types_physics_wrappers

use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper

implicit none

private

    type, public :: Physics_Wrapper
        type(Environment_Wrapper) :: environment
        type(Mixture_Wrapper) :: mixture
        type(Short_Interactions_Wrapper) :: short_interactions
        type(Dipolar_Interactions_Wrapper) :: dipolar_interactions
    end type Physics_Wrapper

end module types_physics_wrappers
