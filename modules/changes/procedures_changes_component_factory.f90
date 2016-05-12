module procedures_changes_component_factory

use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use procedures_moved_coordinates_factory, only: moved_coordinates_create => create, &
    moved_coordinates_destroy => destroy
use module_move_tuning, only: Concrete_Move_Tuning_Parameters
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use procedures_move_tuner_factory, only: move_tuner_create_translation => create_translation, &
    move_tuner_create_rotation => create_rotation, move_tuner_destroy => destroy
use procedures_component_exchange_factory, only: component_exchange_create, &
    component_exchange_destroy
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private
public :: changes_component_create, changes_component_destroy

contains

    subroutine changes_component_create(component, periodic_box, mixture_component, &
        tuning_parameters, tuner_parameters, num_tuning_steps, generating_data, prefix)
        type(Changes_Component_Wrapper), intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: mixture_component
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        type(Concrete_Move_Tuner_Parameters), intent(in) :: tuner_parameters
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        call moved_coordinates_create(component%translated_positions, periodic_box, &
            mixture_component%positions, tuning_parameters, generating_data, prefix)
        call move_tuner_create_translation(component%translation_tuner, component%&
            translated_positions, tuner_parameters, num_tuning_steps)
        call moved_coordinates_create(component%rotated_orientations, mixture_component%&
            orientations, tuning_parameters, generating_data, prefix)
        call move_tuner_create_rotation(component%rotation_tuner, component%&
            rotated_orientations, tuner_parameters, num_tuning_steps)
        call component_exchange_create(component%exchange, mixture_component)
    end subroutine changes_component_create

    subroutine changes_component_destroy(component)
        type(Changes_Component_Wrapper), intent(inout) :: component

        call component_exchange_destroy(component%exchange)
        call move_tuner_destroy(component%rotation_tuner)
        call moved_coordinates_destroy(component%rotated_orientations)
        call move_tuner_destroy(component%translation_tuner)
        call moved_coordinates_destroy(component%translated_positions)
    end subroutine changes_component_destroy

end module procedures_changes_component_factory
