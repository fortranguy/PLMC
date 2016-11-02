module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_input_prefixes, only: changes_prefix
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use procedures_environment_inquirers, only: property_total_volume_can_change => &
    total_volume_can_change
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_inquirers, only: component_can_translate, component_can_rotate, &
    component_can_exchange
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use module_move_tuning, only: Concrete_Move_Tuning_Parameters
use procedures_changed_boxes_size_factory, only: changed_boxes_size_create => create, &
    changed_boxes_size_destroy => destroy
use procedures_exchanged_boxes_size_factory, only: exchanged_boxes_size_create => create, &
    exchanged_boxes_size_destroy => destroy
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use procedures_move_tuner_factory, only: move_tuner_create_boxes_size => create_boxes_size, &
    move_tuner_destroy => destroy
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use procedures_coordinates_copier_factory, only: coordinates_copier_create_position => &
    create_position, coordinates_copier_create_orientation => create_orientation, &
    coordinates_copier_destroy => destroy
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_component_factory, only: changes_component_create => create, &
    changes_component_destroy => destroy
use types_changes_wrapper, only: Changes_Wrapper

implicit none

private
public :: create, destroy, set_can_translate, set_can_rotate, set_can_exchange

interface create
    module procedure :: create_all
    module procedure :: create_components
end interface create

interface destroy
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface destroy

contains

    !> @todo Review the compatibility with GEMC
    subroutine create_all(changes, environment, components, num_tuning_steps, generating_data)
        type(Changes_Wrapper), intent(out) :: changes
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:, :)
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data

        type(Concrete_Move_Tuning_Parameters) :: box_size_tuning_parameters, &
            components_tuning_parameters
        type(Concrete_Move_Tuner_Parameters) :: box_size_tuner_parameters, &
            components_tuner_parameters
        logical, dimension(size(components, 1), size(components, 2)) :: have_positions, &
            have_orientations, can_exchange
        logical :: total_volume_can_change
        logical :: some_boxes_size_can_change, some_components_have_coordinates

        total_volume_can_change = property_total_volume_can_change(environment%beta_pressure)
        some_boxes_size_can_change = total_volume_can_change .or. &
            size(environment%periodic_boxes) > 1
        call set_tuning_parameters(box_size_tuning_parameters, num_tuning_steps, &
            some_boxes_size_can_change, generating_data, changes_prefix//"Boxes.")
        call changed_boxes_size_create(changes%changed_boxes_size, environment%periodic_boxes, &
            box_size_tuning_parameters, total_volume_can_change, generating_data, changes_prefix//&
            "Boxes.")
        call exchanged_boxes_size_create(changes%exchanged_boxes_size, environment%periodic_boxes, &
            box_size_tuning_parameters, generating_data, changes_prefix//"Boxes.")
        call set_tuner_parameters(box_size_tuner_parameters, num_tuning_steps, &
            some_boxes_size_can_change, generating_data, changes_prefix//"Boxes.")
        call move_tuner_create_boxes_size(changes%boxes_size_change_tuner, changes%&
            changed_boxes_size, box_size_tuner_parameters, num_tuning_steps)
        call move_tuner_create_boxes_size(changes%boxes_size_exchange_tuner, changes%&
            exchanged_boxes_size, box_size_tuner_parameters, num_tuning_steps)

        call set_have_positions(have_positions, components)
        call set_have_orientations(have_orientations, components)

        some_components_have_coordinates = any(have_positions) .or. any(have_orientations)
        call set_tuning_parameters(components_tuning_parameters, num_tuning_steps, &
            some_components_have_coordinates, generating_data, changes_prefix//"Components.")
        call set_tuner_parameters(components_tuner_parameters, num_tuning_steps, &
            some_components_have_coordinates, generating_data, changes_prefix//"Components.")

        call create(changes%components, environment%periodic_boxes, components, &
            components_tuning_parameters, components_tuner_parameters, num_tuning_steps, &
            generating_data, changes_prefix//"Components.")

        call set_can_exchange(can_exchange, components)
        call random_coordinates_create(changes%random_positions, environment%accessible_domains, &
            have_positions, can_exchange .or. size(environment%periodic_boxes) > 1)
        call random_coordinates_create(changes%random_orientation, have_orientations, can_exchange)

        call coordinates_copier_create_position(changes%position_copiers, changes%random_positions,&
            have_positions, can_exchange)
        call coordinates_copier_create_orientation(changes%orientation_copier, changes%&
            random_orientation, have_orientations, can_exchange)
    end subroutine create_all

    subroutine create_components(components, periodic_boxes, mixture_components, &
        components_tuning_parameters, components_tuner_parameters, num_tuning_steps, &
        generating_data, prefix)
        type(Changes_Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(Component_Wrapper), intent(in) :: mixture_components(:, :)
        type(Concrete_Move_Tuning_Parameters) :: components_tuning_parameters
        type(Concrete_Move_Tuner_Parameters) :: components_tuner_parameters
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box, i_component
        type(Concrete_Number_to_String) :: string

        allocate(components(size(mixture_components, 1), size(mixture_components, 2)))
        do i_box = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                call changes_component_create(components(i_component, i_box), &
                    periodic_boxes(i_box), mixture_components(i_component, i_box), &
                    components_tuning_parameters, components_tuner_parameters, num_tuning_steps, &
                    generating_data, prefix//"Component "//string%get(i_component)//".")
            end do
        end do
    end subroutine create_components

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
        logical, intent(inout) :: can_translate(:, :)
        type(Changes_Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_translate, 2)
            do i_component = 1, size(can_translate, 1)
                can_translate(i_component, i_box) = &
                    component_can_translate(components(i_component, i_box)%translated_positions)
            end do
        end do
    end subroutine set_can_translate

    subroutine set_can_rotate(can_rotate, components)
        logical, intent(inout) :: can_rotate(:, :)
        type(Changes_Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_rotate, 2)
            do i_component = 1, size(can_rotate, 1)
                can_rotate(i_component, i_box) = &
                    component_can_rotate(components(i_component, i_box)%rotated_orientations)
            end do
        end do
    end subroutine set_can_rotate

    subroutine set_can_exchange(can_exchange, components)
        logical, intent(inout) :: can_exchange(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(can_exchange, 2)
            do i_component = 1, size(can_exchange, 1)
                can_exchange(i_component, i_box) = &
                    component_can_exchange(components(i_component, i_box)%chemical_potential)
            end do
        end do
    end subroutine set_can_exchange

    subroutine destroy_all(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        call coordinates_copier_destroy(changes%orientation_copier)
        call coordinates_copier_destroy(changes%position_copiers)
        call random_coordinates_destroy(changes%random_orientation)
        call random_coordinates_destroy(changes%random_positions)
        call destroy_components(changes%components)
        call move_tuner_destroy(changes%boxes_size_exchange_tuner)
        call move_tuner_destroy(changes%boxes_size_change_tuner)
        call exchanged_boxes_size_destroy(changes%exchanged_boxes_size)
        call changed_boxes_size_destroy(changes%changed_boxes_size)
    end subroutine destroy_all

    subroutine destroy_components(components)
        type(Changes_Component_Wrapper), allocatable, intent(inout) :: components(:, :)

        integer :: i_box
        integer :: i_component

        if (allocated(components)) then
            do i_box = size(components, 2), 1, -1
                do i_component = size(components, 1), 1, -1
                    call changes_component_destroy(components(i_component, i_box))
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

end module procedures_changes_factory
