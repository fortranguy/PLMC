module procedures_physical_model_factory

use json_module, only: json_file
use procedures_environment_factory, only: environment_create => create, environment_destroy => &
    destroy
use procedures_mixture_factory, only: mixture_create => create, mixture_destroy => destroy
use procedures_short_interactions_factory, only: short_interactions_create, &
    short_interactions_destroy
use procedures_dipolar_interactions_factory, only: dipolar_interactions_create, &
    dipolar_interactions_destroy
use procedures_dipolar_interactions_facades_factory, only: dipolar_interactions_facades_create => &
    create, dipolar_interactions_facades_destroy => destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: create_generating, create_exploring, destroy

contains

    subroutine create_generating(physical_model, generating_data)
        type(Physical_Model_Wrapper), intent(out) :: physical_model
        type(json_file), intent(inout) :: generating_data

        call environment_create(physical_model%environment, generating_data)
        call mixture_create(physical_model%mixture, physical_model%environment, generating_data)
        call short_interactions_create(physical_model%short_interactions, physical_model%&
            environment, physical_model%mixture, generating_data)
        call dipolar_interactions_create(physical_model%dipolar_interactions_dynamic, &
            physical_model%gemc_dipolar_interactions_static, physical_model%environment, physical_model%&
            mixture, generating_data)
        call dipolar_interactions_facades_create(physical_model%dipolar_interactions_facades, &
            physical_model%environment, physical_model%mixture%gemc_components, physical_model%&
            dipolar_interactions_dynamic, physical_model%gemc_dipolar_interactions_static)
    end subroutine create_generating

    subroutine create_exploring(physical_model, generating_data, exploring_data, unique_box)
        type(Physical_Model_Wrapper), intent(out) :: physical_model
        type(json_file), intent(inout) :: generating_data, exploring_data
        logical, optional, intent(in) :: unique_box

        call environment_create(physical_model%environment, generating_data, unique_box)
        call mixture_create(physical_model%mixture, physical_model%environment, generating_data)
        call short_interactions_create(physical_model%short_interactions, physical_model%&
            environment, physical_model%mixture, generating_data, exploring_data)
        call dipolar_interactions_create(physical_model%dipolar_interactions_dynamic, &
            physical_model%gemc_dipolar_interactions_static, physical_model%environment, physical_model%&
            mixture, generating_data)
        call dipolar_interactions_facades_create(physical_model%dipolar_interactions_facades, &
            physical_model%environment, physical_model%mixture%gemc_components, physical_model%&
            dipolar_interactions_dynamic, physical_model%gemc_dipolar_interactions_static, &
            exploring_data)
    end subroutine create_exploring

    subroutine destroy(physical_model)
        type(Physical_Model_Wrapper), intent(inout) :: physical_model

        call dipolar_interactions_facades_destroy(physical_model%dipolar_interactions_facades)
        call dipolar_interactions_destroy(physical_model%dipolar_interactions_dynamic, &
            physical_model%gemc_dipolar_interactions_static)
        call short_interactions_destroy(physical_model%short_interactions)
        call mixture_destroy(physical_model%mixture)
        call environment_destroy(physical_model%environment)
    end subroutine destroy

end module procedures_physical_model_factory
