module procedures_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use class_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use types_component_wrapper, only: Component_Wrapper
use module_change_tuning, only: Concrete_Change_Tuning_Parameters
use class_change_tuner, only: Concrete_Change_Tuner_Parameters
use types_changes_wrapper, only: Changes_Wrapper
use procedures_changes_component_factory, only: changes_component_create, changes_component_destroy
use module_plmc_iterations, only: num_tuning_steps

implicit none

private
public :: changes_create, changes_destroy

interface changes_create
    module procedure :: create_all
end interface changes_create

interface changes_destroy
    module procedure :: destroy_all
end interface changes_destroy

contains

    subroutine create_all(changes, periodic_box, components, input_data, prefix)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        type(Concrete_Change_Tuning_Parameters) :: tuning_parameters
        type(Concrete_Change_Tuner_Parameters) :: tuner_parameters
        type(Concrete_Number_to_String) :: string
        integer :: i_component

        call set_tuning_parameters(tuning_parameters, input_data, prefix//"Components.")
        call set_tuner_parameters(tuner_parameters, input_data, prefix//"Components.")
        allocate(changes%components(size(components)))
        do i_component = 1, size(changes%components)
            call changes_component_create(changes%components(i_component), periodic_box, &
                components(i_component), tuning_parameters, tuner_parameters, input_data, &
                prefix//"Component "//string%get(i_component)&
                //".")
        end do
    end subroutine create_all

    subroutine set_tuning_parameters(parameters, input_data, prefix)
        type(Concrete_Change_Tuning_Parameters), intent(out) :: parameters
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (num_tuning_steps > 0) then
            data_field = prefix//"increase factor"
            call input_data%get(data_field, parameters%increase_factor, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"maximum increase factor"
            call input_data%get(data_field, parameters%increase_factor_max, data_found)
            call check_data_found(data_field, data_found)
            deallocate(data_field)
        else
            parameters%increase_factor = 1._DP
            parameters%increase_factor_max = parameters%increase_factor
        end if
    end subroutine set_tuning_parameters

    subroutine set_tuner_parameters(parameters, input_data, prefix)
        type(Concrete_Change_Tuner_Parameters), intent(out) :: parameters
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (num_tuning_steps > 0) then
            data_field = prefix//"accumulation period"
            call input_data%get(data_field, parameters%accumulation_period, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"wanted success ratio"
            call input_data%get(data_field, parameters%wanted_success_ratio, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"tolerance"
            call input_data%get(data_field, parameters%tolerance, data_found)
            call check_data_found(data_field, data_found)
        else
            parameters%accumulation_period = 0
            parameters%wanted_success_ratio = 0._DP
            parameters%tolerance = 0._DP
        end if
    end subroutine set_tuner_parameters

    subroutine destroy_all(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        integer :: i_component

        if (allocated(changes%components)) then
            do i_component = size(changes%components), 1, -1
                call changes_component_destroy(changes%components(i_component))
            end do
            deallocate(changes%components)
        end if
    end subroutine destroy_all

end module procedures_changes_factory
