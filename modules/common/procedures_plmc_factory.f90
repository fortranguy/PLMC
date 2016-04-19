module procedures_plmc_factory

use data_prefixes, only:environment_prefix, mixture_prefix, changes_prefix, &
    short_interactions_prefix, dipolar_interactions_prefix, writers_prefix
use json_module, only: json_file
use procedures_command_arguments, only: create_filename_from_argument
use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_create, environment_destroy
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_factory, only: mixture_create, mixture_destroy, mixture_set
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_factory, only: short_interactions_create, &
    short_interactions_destroy
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipolar_interactions_factory, only: dipolar_interactions_create, &
    dipolar_interactions_destroy
use module_changes_success, only: Concrete_Changes_Success, Concrete_Changes_Counter, &
    Concrete_Switch_Counters, &
    changes_counter_reset, changes_counter_set, switches_counters_reset, switches_counters_set
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: changes_create, changes_destroy
use procedures_metropolis_algorithms_factory, only: metropolis_algorithms_create, &
    metropolis_algorithms_set, metropolis_algorithms_destroy
use types_observables_wrapper, only: Observables_Wrapper
use procedures_observables_factory, only: observables_create, observables_destroy
use types_readers_wrapper, only: Component_Readers_wrapper, Readers_Wrapper
use procedures_readers_factory, only: readers_create, readers_destroy
use types_writers_wrapper, only: Writers_Wrapper
use procedures_writers_factory, only: writers_create, writers_destroy
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper, &
    Metropolis_Algorithm_Pointer
use procedures_property_inquirers, only: use_walls, component_is_dipolar, &
    component_can_change

implicit none

private
public :: plmc_create, plmc_set, plmc_destroy

interface plmc_create
    module procedure :: create_input_data
    module procedure :: create_environment
    module procedure :: create_mixture
    module procedure :: create_short_interactions
    module procedure :: create_dipolar_interactions
    module procedure :: create_changes
    module procedure :: create_metropolis
    module procedure :: create_observables
    module procedure :: create_readers
    module procedure :: create_writers
end interface plmc_create

interface plmc_set
    module procedure :: set_mixture_coordinates
    module procedure :: set_metropolis
    module procedure :: tune_change_components
    module procedure :: set_success_and_reset_counter
end interface plmc_set

interface plmc_destroy
    module procedure :: destroy_writers
    module procedure :: destroy_readers
    module procedure :: destroy_observables
    module procedure :: destroy_metropolis
    module procedure :: destroy_changes
    module procedure :: destroy_dipolar_interactions
    module procedure :: destroy_short_interactions
    module procedure :: destroy_mixture
    module procedure :: destroy_environment
    module procedure :: destroy_input_data
end interface plmc_destroy

contains

    subroutine create_input_data(input_data, i_argument)
        type(json_file), intent(out) :: input_data
        integer, intent(in) :: i_argument

        character(len=:), allocatable :: data_filename

        call create_filename_from_argument(data_filename, i_argument)
        call input_data%load_file(filename = data_filename)
    end subroutine create_input_data

    subroutine destroy_input_data(input_data)
        type(json_file), intent(inout) :: input_data

        call input_data%destroy()
    end subroutine destroy_input_data

    subroutine create_environment(environment, input_data)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: input_data

        call environment_create(environment, input_data, environment_prefix)
    end subroutine create_environment

    subroutine destroy_environment(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call environment_destroy(environment)
    end subroutine destroy_environment

    subroutine create_mixture(mixture, environment, input_data)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(Environment_Wrapper), intent(in) :: environment
        type(json_file), intent(inout) :: input_data

        call mixture_create(mixture, environment, input_data, mixture_prefix)
    end subroutine create_mixture

    subroutine destroy_mixture(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        call mixture_destroy(mixture)
    end subroutine destroy_mixture

    subroutine set_mixture_coordinates(components_readers, input_data)
        type(Component_Readers_wrapper), intent(in) :: components_readers(:)
        type(json_file), intent(inout) :: input_data

        call mixture_set(components_readers, input_data, mixture_prefix)
    end subroutine set_mixture_coordinates

    subroutine create_short_interactions(short_interactions, environment, mixture, input_data)
        type(Short_Interactions_Wrapper), intent(out) :: short_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data

        call short_interactions_create(short_interactions, environment, mixture, input_data, &
            short_interactions_prefix)
    end subroutine create_short_interactions

    subroutine destroy_short_interactions(short_interactions)
        type(Short_Interactions_Wrapper), intent(inout) :: short_interactions

        call short_interactions_destroy(short_interactions)
    end subroutine destroy_short_interactions

    subroutine create_dipolar_interactions(dipolar_interactions, environment, mixture, input_data)
        type(Dipolar_Interactions_Wrapper), intent(out) :: dipolar_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data

        call dipolar_interactions_create(dipolar_interactions, environment, mixture, input_data, &
            dipolar_interactions_prefix)
    end subroutine create_dipolar_interactions

    subroutine destroy_dipolar_interactions(dipolar_interactions)
        type(Dipolar_Interactions_Wrapper), intent(inout) :: dipolar_interactions

        call dipolar_interactions_destroy(dipolar_interactions)
    end subroutine destroy_dipolar_interactions

    subroutine create_changes(changes, periodic_box, components, input_data)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: input_data

        call changes_create(changes, periodic_box, components, input_data, changes_prefix)
    end subroutine create_changes

    subroutine tune_change_components(tuned, i_step, change_components, changes_sucesses)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Component_Wrapper), intent(inout) :: change_components(:)
        type(Concrete_Changes_Success), intent(in) :: changes_sucesses(:)

        logical :: move_tuned(size(change_components)), rotation_tuned(size(change_components))
        integer :: i_component

        do i_component = 1, size(change_components)
            call change_components(i_component)%move_tuner%tune(move_tuned(i_component), i_step, &
                changes_sucesses(i_component)%move)
            call change_components(i_component)%rotation_tuner%tune(rotation_tuned(i_component), &
                i_step, changes_sucesses(i_component)%rotation)
        end do
        tuned = all(move_tuned) .and. all(rotation_tuned)
    end subroutine tune_change_components

    subroutine destroy_changes(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call changes_destroy(changes)
    end subroutine destroy_changes

    subroutine create_metropolis(metropolis, environment, mixture, changes, short_interactions, &
        dipolar_interactions)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Wrapper), intent(in) :: changes
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        call metropolis_algorithms_create(metropolis, environment, mixture, changes%components, &
            short_interactions, dipolar_interactions)
    end subroutine create_metropolis

    subroutine set_metropolis(metropolis, mixture, changes)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Wrapper), intent(in) :: changes

        call metropolis_algorithms_set(metropolis, mixture%components, changes%components)
    end subroutine set_metropolis

    subroutine destroy_metropolis(metropolis)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis

        call metropolis_algorithms_destroy(metropolis)
    end subroutine destroy_metropolis

    subroutine create_observables(observables, components)
        type(Observables_Wrapper), intent(out) :: observables
        type(Component_Wrapper), intent(in) :: components(:)

        call observables_create(observables, components)
    end subroutine create_observables

    subroutine set_success_and_reset_counter(observables)
        type(Observables_Wrapper), intent(inout) :: observables

        call changes_counter_set(observables%changes_sucesses, observables%changes_counters)
        call changes_counter_reset(observables%changes_counters)

        call switches_counters_set(observables%switches_successes, observables%switches_counters)
        call switches_counters_reset(observables%switches_counters)
    end subroutine set_success_and_reset_counter

    subroutine destroy_observables(observables)
        type(Observables_Wrapper), intent(out) :: observables

        call observables_destroy(observables)
    end subroutine destroy_observables

    subroutine create_readers(readers, environment, components)
        type(Readers_Wrapper), intent(out) :: readers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)

        call readers_create(readers, environment%periodic_box, components)
    end subroutine create_readers

    subroutine destroy_readers(readers)
        type(Readers_Wrapper), intent(inout) :: readers

        call readers_destroy(readers)
    end subroutine destroy_readers

    subroutine create_writers(writers, environment, components, changes, short_interactions, &
        dipolar_interactions, input_data)
        type(Writers_Wrapper), intent(out) :: writers
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Wrapper), intent(in) :: changes
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions
        type(json_file), intent(inout) :: input_data

        call writers_create(writers, environment, short_interactions%wall_pairs, components, &
            changes%components, short_interactions%components_pairs, dipolar_interactions%&
            real_pairs, input_data, writers_prefix)
    end subroutine create_writers

    subroutine destroy_writers(writers)
        type(Writers_Wrapper), intent(inout) :: writers

        call writers_destroy(writers)
    end subroutine destroy_writers

end module procedures_plmc_factory
