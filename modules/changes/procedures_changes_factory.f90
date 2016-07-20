module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use types_component_wrapper, only: Component_Wrapper
use procedures_mixture_factory, only: set_have_positions, set_have_orientations
use module_move_tuning, only: Concrete_Move_Tuning_Parameters
use types_move_tuner_parameters, only: Concrete_Move_Tuner_Parameters
use procedures_random_coordinates_factory, only: random_coordinates_create => create, &
    random_coordinates_destroy => destroy
use procedures_coordinates_copier_factory, only: coordinates_copier_create_position => &
    create_position, coordinates_copier_create_orientation => create_orientation, &
    coordinates_copier_destroy => destroy
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use procedures_changes_component_factory, only: changes_component_create, changes_component_destroy
use types_changes_wrapper, only: Changes_Wrapper
use procedures_property_inquirers, only: component_can_translate, component_can_rotate, &
    component_can_exchange

implicit none

private
public :: changes_create, changes_destroy, set_can_translate, set_can_rotate, set_can_exchange

contains

    subroutine changes_create(changes, periodic_box, accessible_domain, components, &
        num_tuning_steps, generating_data, prefix)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domain
        type(Component_Wrapper), intent(in) :: components(:)
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        type(Concrete_Move_Tuning_Parameters) :: tuning_parameters
        type(Concrete_Move_Tuner_Parameters) :: tuner_parameters
        logical :: can_exchange(size(components))
        logical :: have_positions(size(components)), have_orientations(size(can_exchange))
        type(Concrete_Number_to_String) :: string
        integer :: i_component

        call set_tuning_parameters(tuning_parameters, num_tuning_steps, generating_data, &
            prefix//"Components.")
        call set_tuner_parameters(tuner_parameters, num_tuning_steps, generating_data, &
            prefix//"Components.")
        allocate(changes%components(size(components)))
        do i_component = 1, size(changes%components)
            call changes_component_create(changes%components(i_component), periodic_box, &
                components(i_component), tuning_parameters, tuner_parameters, num_tuning_steps, &
                generating_data, prefix//"Component "//string%get(i_component)//".")
        end do

        call set_can_exchange(can_exchange, components)
        call set_have_positions(have_positions, components)
        call random_coordinates_create(changes%random_position, accessible_domain, can_exchange, &
        have_positions)
        call set_have_orientations(have_orientations, components)
        call random_coordinates_create(changes%random_orientation, can_exchange, have_orientations)
        call coordinates_copier_create_position(changes%position_copier, changes%random_position, &
            have_positions)
        call coordinates_copier_create_orientation(changes%orientation_copier, changes%&
            random_orientation, have_orientations)
    end subroutine changes_create

    subroutine set_tuning_parameters(parameters, num_tuning_steps, generating_data, prefix)
        type(Concrete_Move_Tuning_Parameters), intent(out) :: parameters
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (num_tuning_steps > 0) then
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

    subroutine set_tuner_parameters(parameters, num_tuning_steps, generating_data, prefix)
        type(Concrete_Move_Tuner_Parameters), intent(out) :: parameters
        integer, intent(in) :: num_tuning_steps
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (num_tuning_steps > 0) then
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
    end subroutine changes_destroy

end module procedures_changes_factory
