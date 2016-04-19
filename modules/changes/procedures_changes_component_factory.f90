module procedures_changes_component_factory

use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use procedures_changed_coordinates_factory, only: changed_coordinates_create, &
    changed_coordinates_destroy
use module_change_tuning, only: Concrete_Change_Tuning_Parameters
use types_change_tuner_parameters, only: Concrete_Change_Tuner_Parameters
use procedures_change_tuner_factory, only: change_tuner_create_move, change_tuner_create_rotation, &
    change_tuner_destroy
use procedures_component_exchange_factory, only: component_exchange_create, &
    component_exchange_destroy
use types_changes_component_wrapper, only: Changes_Component_Wrapper

implicit none

private
public :: changes_component_create, changes_component_destroy

contains

    subroutine changes_component_create(component, periodic_box, mixture_component, &
        tuning_parameters, tuner_parameters, input_data, prefix)
        type(Changes_Component_Wrapper), intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: mixture_component
        type(Concrete_Change_Tuning_Parameters), intent(in) :: tuning_parameters
        type(Concrete_Change_Tuner_Parameters), intent(in) :: tuner_parameters
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call changed_coordinates_create(component%moved_positions, periodic_box, mixture_component%&
            positions, tuning_parameters, input_data, prefix)
        call change_tuner_create_move(component%move_tuner, component%moved_positions, &
            tuner_parameters)
        call changed_coordinates_create(component%rotated_orientations, mixture_component%&
            orientations, tuning_parameters, input_data, prefix)
        call change_tuner_create_rotation(component%rotation_tuner, component%&
            rotated_orientations, tuner_parameters)
        call component_exchange_create(component%exchange, mixture_component)
    end subroutine changes_component_create

    subroutine changes_component_destroy(component)
        type(Changes_Component_Wrapper), intent(inout) :: component

        call component_exchange_destroy(component%exchange)
        call change_tuner_destroy(component%rotation_tuner)
        call changed_coordinates_destroy(component%rotated_orientations)
        call change_tuner_destroy(component%move_tuner)
        call changed_coordinates_destroy(component%moved_positions)
    end subroutine changes_component_destroy

end module procedures_changes_component_factory
