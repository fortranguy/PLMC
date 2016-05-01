program plmc_widom

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use types_observables_wrapper, only: Observables_Wrapper
use types_readers_wrapper, only: Readers_Wrapper
use types_writers_wrapper, only: Writers_Wrapper
use procedures_plmc_factory, only: plmc_create, plmc_destroy, plmc_set
use procedures_plmc_reset, only: plmc_reset
use procedures_plmc_help, only: plmc_widom_catch_help

implicit none

    type(Environment_Wrapper) :: environment
    type(Mixture_Wrapper) :: mixture
    type(Changes_Wrapper) :: changes
    type(Short_Interactions_Wrapper) :: short_interactions
    type(Dipolar_Interactions_Wrapper) :: dipolar_interactions
    type(Metropolis_Algorithms_Wrapper) :: metropolis_algorithms
    type(Observables_Wrapper) :: observables
    type(Readers_Wrapper) :: readers
    type(Writers_Wrapper) :: writers

    type(json_file) :: input_data, post_data
    integer :: i_step

    call plmc_widom_catch_help()
    call plmc_create(input_data, command_argument_count() - 1)
    !call plmc_create(environment, mixture, short_interactions, dipolar_interactions, changes, &
    !    metropolis_algorithms, observables, readers, writers, input_data)
    call plmc_set(readers%components, input_data)
    call plmc_destroy(input_data)

    call plmc_create(post_data, command_argument_count())
    call plmc_destroy(post_data)

    !call plmc_destroy(environment, mixture, short_interactions, dipolar_interactions, changes, &
    !    metropolis_algorithms, observables, readers, writers)

end program plmc_widom
