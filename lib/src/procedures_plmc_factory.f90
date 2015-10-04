module procedures_plmc_factory

use data_constants, only: num_components
use json_module, only: json_file, json_initialize
use procedures_command_arguments, only: set_filename_from_argument
use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use types_particles_wrapper, only: Mixture_Wrapper, Particles_Wrapper
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_factory, only: changes_factory_create, changes_factory_destroy
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use types_observable_writers_wrapper, only: Mixture_Observable_Writers_Wrapper
use procedures_observable_writers_factory, only: observable_writers_factory_create, &
    observable_writers_factory_destroy

implicit none

private
public :: plmc_load, plmc_create, plmc_destroy

interface plmc_load
    module procedure :: load_input_data
end interface plmc_load

interface plmc_create
    module procedure :: create_environment
    module procedure :: create_mixture
    module procedure :: create_changes
    module procedure :: create_short_potentials
    module procedure :: create_observable_writers
end interface plmc_create

interface plmc_destroy
    module procedure :: destroy_observable_writers
    module procedure :: destroy_short_potentials
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

        call environment_factory_create(environment, input_data, "Environment.")
    end subroutine create_environment

    subroutine destroy_environment(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call environment_factory_destroy(environment)
    end subroutine destroy_environment

    subroutine create_mixture(mixture, input_data, environment)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(json_file), intent(inout) :: input_data
        type(Environment_Wrapper), intent(in) :: environment

        call particles_factory_create(mixture%components(1), input_data, "Mixture.Component 1.", &
            environment)
        call particles_factory_create(mixture%components(2), input_data, "Mixture.Component 2.", &
            environment)
        call particles_factory_create(mixture%inter_diameter, mixture%components(1)%diameter, &
            mixture%components(2)%diameter, input_data, "Mixture.Inter 12.")
    end subroutine create_mixture

    subroutine destroy_mixture(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        call particles_factory_destroy(mixture%inter_diameter)
        call particles_factory_destroy(mixture%components(2))
        call particles_factory_destroy(mixture%components(1))
    end subroutine destroy_mixture

    subroutine create_changes(changes, input_data, periodic_box, components)
        type(Changes_Wrapper), intent(out) :: changes(num_components)
        type(json_file), intent(inout) :: input_data
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Particles_Wrapper), intent(in) :: components(num_components)

        call changes_factory_create(changes(1), input_data, "Changes.Component 1.", periodic_box, &
            components(1))
        call changes_factory_create(changes(2), input_data, "Changes.Component 2.", periodic_box, &
            components(2))
    end subroutine create_changes

    subroutine destroy_changes(changes)
        type(Changes_Wrapper), intent(inout) :: changes(2)

        call changes_factory_destroy(changes(1))
        call changes_factory_destroy(changes(2))
    end subroutine destroy_changes

    subroutine create_short_potentials(short_potentials, input_data, environment, mixture)
        type(Mixture_Short_Potentials_Wrapper), intent(out) :: short_potentials
        type(json_file), intent(inout) :: input_data
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture

        call short_potential_factory_create(short_potentials%intras(1), input_data, &
            "Short Potentials.Component 1.", environment, mixture%components(1))
        call short_potential_factory_create(short_potentials%intras(2), input_data, &
            "Short Potentials.Component 2.", environment, mixture%components(2))
        call short_potential_factory_create(short_potentials%inter_micro, input_data, &
            "Short Potentials.Inter 12.", mixture%inter_diameter)
        call short_potential_factory_create(short_potentials%inters(1), &
            short_potentials%inter_micro, input_data, "Short Potentials.Inter 12.", &
            environment%periodic_box, mixture%components(1)%positions)
        call short_potential_factory_create(short_potentials%inters(2), &
            short_potentials%inter_micro, input_data, "Short Potentials.Inter 12.", &
            environment%periodic_box, mixture%components(2)%positions)
    end subroutine create_short_potentials

    subroutine destroy_short_potentials(short_potentials)
        type(Mixture_Short_Potentials_Wrapper), intent(inout) :: short_potentials

        call short_potential_factory_destroy(short_potentials%inters(2))
        call short_potential_factory_destroy(short_potentials%inters(1))
        call short_potential_factory_destroy(short_potentials%inter_micro)
        call short_potential_factory_destroy(short_potentials%intras(2))
        call short_potential_factory_destroy(short_potentials%intras(1))
    end subroutine destroy_short_potentials

    subroutine create_observable_writers(observable_writers, input_data, mixture, changes)
        type(Mixture_Observable_Writers_Wrapper), intent(out) :: observable_writers
        type(json_file), intent(inout) :: input_data
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        call observable_writers_factory_create(observable_writers%intras(1)%energy, &
            mixture%components(1)%number, mixture%components(1)%wall_diameter, &
            "component_1_energy.out")
        call observable_writers_factory_create(observable_writers%intras(2)%energy, &
            mixture%components(2)%number, mixture%components(2)%wall_diameter, &
            "component_2_energy.out")
        call observable_writers_factory_create(observable_writers%inter_energy, &
            mixture%inter_diameter, "inter_12_energy.out")
        call observable_writers_factory_create(observable_writers%intras(1)%changes, &
            changes(1)%moved_positions, changes(1)%rotated_orientations, &
            changes(1)%particles_exchange, "component_1_success.out")
        call observable_writers_factory_create(observable_writers%intras(2)%changes, &
            changes(2)%moved_positions, changes(2)%rotated_orientations, &
            changes(2)%particles_exchange, "component_2_success.out")
        call observable_writers_factory_create(observable_writers%intras(1)%coordinates, &
            "component_1_coordinates", mixture%components(1)%positions, &
            mixture%components(1)%orientations, input_data, "Monte Carlo.")
        call observable_writers_factory_create(observable_writers%intras(2)%coordinates, &
            "component_2_coordinates", mixture%components(2)%positions, &
            mixture%components(2)%orientations, input_data, "Monte Carlo.")
    end subroutine create_observable_writers

    subroutine destroy_observable_writers(observable_writers)
        type(Mixture_Observable_Writers_Wrapper), intent(inout) :: observable_writers
        call observable_writers_factory_destroy(observable_writers%intras(2)%coordinates)
        call observable_writers_factory_destroy(observable_writers%intras(1)%coordinates)
        call observable_writers_factory_destroy(observable_writers%intras(2)%changes)
        call observable_writers_factory_destroy(observable_writers%intras(1)%changes)
        call observable_writers_factory_destroy(observable_writers%inter_energy)
        call observable_writers_factory_destroy(observable_writers%intras(2)%energy)
        call observable_writers_factory_destroy(observable_writers%intras(1)%energy)
    end subroutine destroy_observable_writers

end module procedures_plmc_factory
