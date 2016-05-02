module procedures_physical_model_factory

use data_prefixes, only: environment_prefix, mixture_prefix, short_interactions_prefix, &
    dipolar_interactions_prefix
use json_module, only: json_file
use procedures_environment_factory, only: environment_create, environment_destroy
use procedures_mixture_factory, only: mixture_create, mixture_destroy
use procedures_short_interactions_factory, only: short_interactions_create, &
    short_interactions_destroy
use procedures_dipolar_interactions_factory, only: dipolar_interactions_create, &
    dipolar_interactions_destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: create, destroy

contains

    subroutine create(physical_model, input_data)
        type(Physical_Model_Wrapper), intent(out) :: physical_model
        type(json_file), intent(inout) :: input_data

        call environment_create(physical_model%environment, input_data, environment_prefix)
        call mixture_create(physical_model%mixture, physical_model%environment, input_data, &
            mixture_prefix)
        call short_interactions_create(physical_model%short_interactions, physical_model%&
            environment, physical_model%mixture, input_data, short_interactions_prefix)
        call dipolar_interactions_create(physical_model%dipolar_interactions, physical_model%&
            environment, physical_model%mixture, input_data, dipolar_interactions_prefix)
    end subroutine create

    subroutine destroy(physical_model)
        type(Physical_Model_Wrapper), intent(inout) :: physical_model

        call dipolar_interactions_destroy(physical_model%dipolar_interactions)
        call short_interactions_destroy(physical_model%short_interactions)
        call mixture_destroy(physical_model%mixture)
        call environment_destroy(physical_model%environment)
    end subroutine destroy

end module procedures_physical_model_factory
