module procedures_short_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_potential_expression, only: Abstract_Potential_Expression, &
    Null_Potential_Expression, Lennard_Jones_Expression
use class_pair_potential, only: Abstract_Pair_Potential, &
    Concrete_Pair_Potential, Null_Pair_Potential, Hard_Pair_Potential
use types_short_potential, only: Short_Potential_Wrapper

implicit none

private
public :: allocate_and_set_potential_expression

contains

    subroutine potential_factory_construct(short_potential, input_data, prefix)
        type(Short_Potential_Wrapper), intent(out) :: short_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

    end subroutine potential_factory_construct

    subroutine allocate_and_set_potential_expression(potential_expression, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_potential_expression(potential_expression, input_data, prefix)
        call set_potential_expression(potential_expression, input_data, prefix)
    end subroutine allocate_and_set_potential_expression

    subroutine allocate_potential_expression(potential_expression, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        !if particles exist...
        data_field = prefix//".Potential.name"
        call input_data%get(data_field, potential_name, data_found)
        call test_data_found(data_field, data_found)
        select case(potential_name)
            case ("HS")
                allocate(Null_Potential_Expression :: potential_expression)

            case ("LJ")
                allocate(Lennard_Jones_Expression :: potential_expression)
            case default
                call error_exit(potential_name//" unknown potential_name."//&
                    "Choose between: 'HS' and LJ.")
        end select
        deallocate(potential_name)
        deallocate(data_field)
    end subroutine allocate_potential_expression

    subroutine set_potential_expression(potential_expression, input_data, prefix)
        class(Abstract_Potential_Expression), intent(inout) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: LJ_epsilon, LJ_sigma

        select type(potential_expression)
            type is (Null_Potential_Expression)
                call potential_expression%set()
            type is (Lennard_Jones_Expression)
                data_field = "Potential.epsilon"
                call input_data%get(data_field, LJ_epsilon, data_found)
                call test_data_found(data_field, data_found)
                data_field = "Potential.sigma"
                call input_data%get(data_field, LJ_sigma, data_found)
                call test_data_found(data_field, data_found)
                call potential_expression%set(LJ_epsilon, LJ_sigma)
        end select
        if (allocated(data_field)) deallocate(data_field)
    end subroutine set_potential_expression

    subroutine potential_factory_destroy(short_potential)
        type(Short_Potential_Wrapper), intent(inout) :: short_potential
    end subroutine potential_factory_destroy

end module procedures_short_potential_factory
