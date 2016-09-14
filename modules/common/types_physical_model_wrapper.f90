module types_physical_model_wrapper

use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper

implicit none

private

    type, public :: Physical_Model_Wrapper
        type(Environment_Wrapper) :: environment
        type(Mixture_Wrapper) :: mixture
        type(Short_Interactions_Wrapper) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper) :: dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper) :: dipolar_interactions_static
    end type Physical_Model_Wrapper

end module types_physical_model_wrapper
