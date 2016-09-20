module procedures_dipolar_interactions_facade_factory

use json_module, only: json_file
use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: periodicity_is_xyz, box_size_can_change
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_total_moment_factory, only: set_are_dipolar
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use classes_dipolar_interactions_facade, only: Abstract_Dipolar_Interactions_Facade, &
    Scalable_Dipolar_Interactions_Facade, Unscalable_Dipolar_Interactions_Facade, &
    Null_Dipolar_Interactions_Facade
use procedures_exploration_inquirers, only: property_measure_pressure => measure_pressure

implicit none

private
public :: create, destroy

contains

    subroutine create(facade, environment, components, dipolar_interactions_dynamic, &
        dipolar_interactions_static, exploring_data, volume_change_prefix)
        class(Abstract_Dipolar_Interactions_Facade), allocatable, intent(out) :: facade
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        type(Dipolar_Interactions_Dynamic_Wrapper), intent(in) :: dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), intent(in) :: dipolar_interactions_static
        type(json_file), optional, intent(inout) :: exploring_data
        character(len=*), optional, intent(in) :: volume_change_prefix

        logical :: are_dipolar(size(components)), measure_pressure

        call set_are_dipolar(are_dipolar, components)
        if (present(exploring_data) .and. present(volume_change_prefix)) then
            measure_pressure = property_measure_pressure(exploring_data, volume_change_prefix)
        else
            measure_pressure = .false.
        end if
        if ((box_size_can_change(environment%beta_pressure) .or. measure_pressure) .and. &
            any(are_dipolar)) then
            if (periodicity_is_xyz(environment%periodic_box)) then
                allocate(Scalable_Dipolar_Interactions_Facade :: facade)
            else
                allocate(Unscalable_Dipolar_Interactions_Facade :: facade)
            end if
        else
            allocate(Null_Dipolar_Interactions_Facade :: facade)
        end if

        select type (facade)
            type is (Scalable_Dipolar_Interactions_Facade)
                call facade%construct(dipolar_interactions_dynamic, dipolar_interactions_static)
            type is (Unscalable_Dipolar_Interactions_Facade)
                call facade%construct(environment%periodic_box, components, &
                    dipolar_interactions_dynamic, dipolar_interactions_static)
            type is (Null_Dipolar_Interactions_Facade)
            class default
                call error_exit("procedures_dipolar_interactions_facade_factory: create: "//&
                    "facade type is unknown.")
        end select
    end subroutine create

    subroutine destroy(facade)
        class(Abstract_Dipolar_Interactions_Facade), allocatable, intent(inout) :: facade

        if (allocated(facade)) then
            call facade%destroy()
            deallocate(facade)
        end if
    end subroutine destroy

end module procedures_dipolar_interactions_facade_factory
