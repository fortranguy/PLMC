module procedures_plmc_factory

use data_wrappers_prefix, only:environment_prefix, mixture_prefix, changes_prefix, &
    short_interactions_prefix, long_interactions_prefix
use json_module, only: json_file, json_initialize
use procedures_command_arguments, only: set_filename_from_argument
use class_periodic_box, only: Abstract_Periodic_Box
use class_walls_potential, only: Abstract_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_create, environment_destroy
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_component_factory, only: component_create, component_destroy !to delete
use procedures_mixture_factory, only: mixture_create, mixture_destroy
use module_changes_success, only: Concrete_Changes_Success, Concrete_Changes_Counter
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: changes_create, changes_destroy
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_short_interactions_factory, only: short_interactions_create, &
    short_interactions_destroy
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use procedures_long_interactions_factory, only: long_interactions_create, long_interactions_destroy
use module_changes_success, only: reset_counter => Concrete_Changes_Counter_reset, &
    set_success => Concrete_Changes_Counter_set
use types_observables_wrapper, only: Observables_Wrapper
use types_writers_wrapper, only: Writers_Wrapper
use procedures_writers_factory, only: writers_create, &
    writers_destroy
use types_metropolis_wrapper, only: Metropolis_Wrapper
use procedures_metropolis_factory, only: metropolis_create, metropolis_set, &
    metropolis_destroy
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use procedures_property_inquirers, only: use_walls, component_is_dipolar, &
    component_can_change

implicit none

private
public :: plmc_load, plmc_create, plmc_set, plmc_destroy

interface plmc_load
    module procedure :: load_input_data
end interface plmc_load

interface plmc_create
    module procedure :: create_environment
    module procedure :: create_mixture
    module procedure :: create_changes
    module procedure :: create_short_interactions
    module procedure :: create_long_interactions
    module procedure :: create_writers
    module procedure :: create_metropolis
end interface plmc_create

interface plmc_set
    module procedure :: set_metropolis
    module procedure :: tune_changes
    module procedure :: set_success_and_reset_counter
end interface plmc_set

interface plmc_destroy
    module procedure :: destroy_metropolis
    module procedure :: destroy_writers
    module procedure :: destroy_long_interactions
    module procedure :: destroy_short_interactions
    module procedure :: destroy_changes
    module procedure :: destroy_mixture
    module procedure :: destroy_environment
end interface plmc_destroy

contains

    subroutine load_input_data(input_data)
        type(json_file), intent(out) :: input_data

        character(len=:), allocatable :: data_filename

        call json_initialize()
        call set_filename_from_argument(data_filename)
        call input_data%load_file(filename = data_filename)
        if (allocated(data_filename)) deallocate(data_filename)
    end subroutine load_input_data

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

    subroutine create_changes(changes, periodic_box, components, input_data)
        type(Changes_Wrapper), intent(out) :: changes(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: input_data

        call changes_create(changes(1), periodic_box, components(1), input_data, &
            changes_prefix//"Component 1.")
        call changes_create(changes(2), periodic_box, components(2), input_data, &
            changes_prefix//"Component 2.")
    end subroutine create_changes

    subroutine destroy_changes(changes)
        type(Changes_Wrapper), intent(inout) :: changes(2)

        call changes_destroy(changes(1))
        call changes_destroy(changes(2))
    end subroutine destroy_changes

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

    subroutine create_long_interactions(long_interactions, environment, mixture, input_data)
        type(Long_Interactions_Wrapper), intent(out) :: long_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data

        call long_interactions_create(long_interactions, environment, mixture, input_data, &
            long_interactions_prefix)
    end subroutine create_long_interactions

    subroutine destroy_long_interactions(long_interactions)
        type(Long_Interactions_Wrapper), intent(inout) :: long_interactions

        call long_interactions_destroy(long_interactions)
    end subroutine destroy_long_interactions

    subroutine create_writers(writers, walls_potential, mixture, changes, &
        input_data)
        type(Writers_Wrapper), intent(out) :: writers
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Wrapper), intent(in) :: changes(:)
        type(json_file), intent(inout) :: input_data

        logical :: wall_used
        logical :: mixture_exists, mixture_is_dipolar
        logical :: component_1_exists, component_1_is_dipolar
        logical :: component_2_exists, component_2_is_dipolar

        wall_used = use_walls(walls_potential)
        component_1_exists = .false.
        !component_1_is_dipolar = component_is_dipolar(mixture%components(1)%dipolar_moments)
        !call writers_create(writers%intras(1)%energy, wall_used, &
        !    component_1_exists, component_1_is_dipolar, "component_1_energy.out")
        !component_2_exists = .false.
        !component_2_is_dipolar = component_is_dipolar(mixture%components(2)%dipolar_moments)
        !call writers_create(writers%intras(2)%energy, wall_used, &
        !    component_2_exists, component_2_is_dipolar, "component_2_energy.out")
        mixture_exists = component_1_exists .and. component_2_exists
        mixture_is_dipolar = component_1_is_dipolar .and. component_2_is_dipolar
        !call writers_create(writers%components_energy, mixture_exists, &
        !    mixture_is_dipolar, "components_12_energy.out")

        !call writers_create(writers%intras(1)%changes, &
        !    changes(1)%moved_positions, changes(1)%rotated_orientations, &
        !    changes(1)%component_exchange, "component_1_success.out")
        !call writers_create(writers%intras(2)%changes, &
        !    changes(2)%moved_positions, changes(2)%rotated_orientations, &
        !    changes(2)%component_exchange, "component_2_success.out")

        !call writers_create(writers%intras(1)%coordinates, &
        !    "component_1_coordinates", mixture%components(1)%positions, &
        !    mixture%components(1)%orientations, input_data, "Monte Carlo.")
        !call writers_create(writers%intras(2)%coordinates, &
        !    "component_2_coordinates", mixture%components(2)%positions, &
        !    mixture%components(2)%orientations, input_data, "Monte Carlo.")
    end subroutine create_writers

    subroutine destroy_writers(writers)
        type(Writers_Wrapper), intent(inout) :: writers

    end subroutine destroy_writers

    subroutine create_metropolis(metropolis, environment, changes)
        type(Metropolis_Wrapper), intent(out) :: metropolis
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(:)

        call metropolis_create(metropolis, environment, changes)
    end subroutine create_metropolis

    subroutine destroy_metropolis(metropolis)
        type(Metropolis_Wrapper), intent(inout) :: metropolis

        call metropolis_destroy(metropolis)
    end subroutine destroy_metropolis

    subroutine set_metropolis(metropolis, components, short_interactions, long_interactions)
        type(Metropolis_Wrapper), intent(inout) :: metropolis
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        call metropolis_set(metropolis, components, short_interactions, long_interactions)
    end subroutine set_metropolis

    subroutine tune_changes(tuned, i_step, changes, changes_sucesses)
        logical, intent(out) :: tuned
        integer, intent(in) :: i_step
        type(Changes_Wrapper), intent(inout) :: changes(:)
        type(Concrete_Changes_Success), intent(in) :: changes_sucesses(:)

        logical :: move_tuned(size(changes)), rotation_tuned(size(changes))
        integer :: i_component

        do i_component = 1, size(changes)
            call changes(i_component)%move_tuner%tune(move_tuned(i_component), i_step, &
                changes_sucesses(i_component)%move)
            call changes(i_component)%rotation_tuner%tune(rotation_tuned(i_component), i_step, &
                changes_sucesses(i_component)%rotation)
        end do
        tuned = all(move_tuned) .and. all(rotation_tuned)
    end subroutine tune_changes

    subroutine set_success_and_reset_counter(changes_sucesses, changes_counters)
        type(Concrete_Changes_Success), intent(out) :: changes_sucesses(:)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters(:)

        call set_success(changes_sucesses, changes_counters)
        call reset_counter(changes_counters)
    end subroutine set_success_and_reset_counter

end module procedures_plmc_factory
