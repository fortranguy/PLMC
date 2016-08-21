module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_input_prefixes, only: environment_prefix
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use module_move_tuning, only: Concrete_Move_Tuning_Parameters
use procedures_changed_box_size_factory, only: changed_box_size_create => create, &
    changed_box_size_destroy => destroy
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use procedures_move_tuner_factory, only: move_tuner_create_box_size_change => &
    create_box_size_change, move_tuner_destroy => destroy
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use procedures_coordinates_copier_factory, only: coordinates_copier_create_position => &
    create_position, coordinates_copier_create_orientation => create_orientation, &
    coordinates_copier_destroy => destroy
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_component_factory, only: changes_component_create, changes_component_destroy
use types_changes_wrapper, only: Changes_Wrapper
use procedures_property_inquirers, only: box_size_can_change, component_can_translate, &
    component_can_rotate, component_can_exchange

implicit none

private
public :: changes_create, changes_destroy, set_can_translate, set_can_rotate, set_can_exchange

contains

    subroutine changes_create(changes, environment, components, num_tuning_steps, generating_data, &
        prefix)
        type(Changes_Wrapper), intent(out) :: changes
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        type(Concrete_Move_Tuning_Parameters) :: box_size_tuning_parameters, &
            components_tuning_parameters
        type(Concrete_Move_Tuner_Parameters) :: box_size_tuner_parameters, &
            components_tuner_parameters
        logical :: can_exchange(size(components))
        logical :: have_positions(size(components)), have_orientations(size(can_exchange))
        logical :: volume_can_change, some_components_have_coordinates
        type(Concrete_Number_to_String) :: string
        integer :: i_component

        volume_can_change = box_size_can_change(environment%beta_pressure)
        call set_tuning_parameters(box_size_tuning_parameters, num_tuning_steps, volume_can_change,&
            generating_data, prefix//"Box Size.")
        call changed_box_size_create(changes%changed_box_size, environment%periodic_box, &
            box_size_tuning_parameters, volume_can_change, generating_data, prefix//"Box Size.")
        call set_tuner_parameters(box_size_tuner_parameters, num_tuning_steps, volume_can_change, &
            generating_data, prefix//"Box Size.")
        call move_tuner_create_box_size_change(changes%box_size_change_tuner, changes%&
            changed_box_size, box_size_tuner_parameters, num_tuning_steps)

        call set_have_positions(have_positions, components)
        call set_have_orientations(have_orientations, components)
        some_components_have_coordinates = any(have_positions) .or. any(have_orientations)
        call set_tuning_parameters(components_tuning_parameters, num_tuning_steps, &
            some_components_have_coordinates, generating_data, prefix//"Components.")
        call set_tuner_parameters(components_tuner_parameters, num_tuning_steps, &
            some_components_have_coordinates, generating_data, prefix//"Components.")
        allocate(changes%components(size(components)))
        do i_component = 1, size(changes%components)
            call changes_component_create(changes%components(i_component), environment%&
                periodic_box, components(i_component), components_tuning_parameters, &
                components_tuner_parameters, num_tuning_steps, generating_data, prefix//&
                "Component "//string%get(i_component)//".")
        end do

        call set_can_exchange(can_exchange, components)
        call random_coordinates_create(changes%random_position, environment%accessible_domain, &
            can_exchange, have_positions)
        call random_coordinates_create(changes%random_orientation, can_exchange, have_orientations)
        call coordinates_copier_create_position(changes%position_copier, changes%random_position, &
            have_positions)
        call coordinates_copier_create_orientation(changes%orientation_copier, changes%&
            random_orientation, have_orientations)
    end subroutine changes_create

    subroutine set_tuning_parameters(parameters, num_tuning_steps, needed, generating_data, prefix)
        type(Concrete_Move_Tuning_Parameters), intent(out) :: parameters
        integer, intent(in) :: num_tuning_steps
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (needed .and. num_tuning_steps > 0) then
            data_field = prefix//"increase factor"
            call generating_data%get(data_field, parameters%increase_factor, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"maximum increase factor"
            call generating_data%get(data_field, parameters%increase_factor_max, data_found)
            call check_data_found(data_field, data_found)
        else
            parameters%increase_factor = 1._DP
            parameters%increase_factor_max = parameters%increase_factor
        end if
    end subroutine set_tuning_parameters

    subroutine set_tuner_parameters(parameters, num_tuning_steps, needed, generating_data, prefix)
        type(Concrete_Move_Tuner_Parameters), intent(out) :: parameters
        integer, intent(in) :: num_tuning_steps
        logical, intent(in) :: needed
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (needed .and. num_tuning_steps > 0) then
            data_field = prefix//"accumulation period"
            call generating_data%get(data_field, parameters%accumulation_period, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"wanted success ratio"
            call generating_data%get(data_field, parameters%wanted_success_ratio, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"tolerance"
            call generating_data%get(data_field, parameters%tolerance, data_found)
            call check_data_found(data_field, data_found)
        else
            parameters%accumulation_period = 0
            parameters%wanted_success_ratio = 0._DP
            parameters%tolerance = 0._DP
        end if
    end subroutine set_tuner_parameters

    subroutine set_can_translate(can_translate, components)
        logical, intent(out) :: can_translate(:)
        type(Changes_Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(can_translate)
            can_translate(i_component) = component_can_translate(components(i_component)%&
                translated_positions)
        end do
    end subroutine set_can_translate

    subroutine set_can_rotate(can_rotate, components)
        logical, intent(out) :: can_rotate(:)
        type(Changes_Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(can_rotate)
            can_rotate(i_component) = component_can_rotate(components(i_component)%&
                rotated_orientations)
        end do
    end subroutine set_can_rotate

    subroutine set_can_exchange(can_exchange, components)
        logical, intent(out) :: can_exchange(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(can_exchange)
            can_exchange(i_component) = component_can_exchange(components(i_component)%&
                chemical_potential)
        end do
    end subroutine set_can_exchange

    subroutine changes_destroy(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        integer :: i_component

        call coordinates_copier_destroy(changes%orientation_copier)
        call coordinates_copier_destroy(changes%position_copier)
        call random_coordinates_destroy(changes%random_orientation)
        call random_coordinates_destroy(changes%random_position)
        if (allocated(changes%components)) then
            do i_component = size(changes%components), 1, -1
                call changes_component_destroy(changes%components(i_component))
            end do
            deallocate(changes%components)
        end if
        call move_tuner_destroy(changes%box_size_change_tuner)
        call changed_box_size_destroy(changes%changed_box_size)
    end subroutine changes_destroy

end module procedures_changes_factory
