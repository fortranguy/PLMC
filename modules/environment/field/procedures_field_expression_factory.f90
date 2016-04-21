module procedures_field_expression_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_permittivity, only: Abstract_Permittivity
use classes_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression, &
    Centered_Plates_Expression, Null_Field_Expression

implicit none

private
public :: create, destroy

contains

    subroutine create(field_expression, permittivity, field_applied, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        class(Abstract_Permittivity), intent(in) :: permittivity
        logical, intent(in) :: field_applied
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_field_expression(field_expression, field_applied, input_data, prefix)
        call set_field_expression(field_expression, permittivity, input_data, prefix)
    end subroutine create

    subroutine allocate_field_expression(field_expression, field_applied, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        logical, intent(in) :: field_applied
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: field_name, data_field
        logical :: data_found

        if (field_applied) then
            data_field = prefix//"External Field.name"
            call input_data%get(data_field, field_name, data_found)
            call check_data_found(data_field, data_found)
            select case (field_name)
                case ("constant")
                    allocate(Constant_Field_Expression :: field_expression)
                case ("plates")
                    allocate(Centered_Plates_Expression :: field_expression)
                case default
                    call error_exit(field_name//" field_name unknown. Choose:"//&
                        " 'constant' or 'plates'.")
            end select
        else
            allocate(Null_Field_Expression :: field_expression)
        end if
    end subroutine allocate_field_expression

    subroutine set_field_expression(field_expression, permittivity, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP), allocatable :: field_vector(:)
        real(DP) :: gap, size_x, surface_density
        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (field_expression)
            type is (Null_Field_Expression)
                call field_expression%set()
            type is (Constant_Field_Expression)
                data_field = prefix//"External Field.vector"
                call input_data%get(data_field, field_vector, data_found)
                call check_data_found(data_field, data_found)
                call field_expression%set(field_vector)
                deallocate(field_vector)
            type is (Centered_Plates_Expression)
                data_field = prefix//"External Field.gap"
                call input_data%get(data_field, gap, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"External Field.size x"
                call input_data%get(data_field, size_x, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"External Field.surface density"
                call input_data%get(data_field, surface_density, data_found)
                call check_data_found(data_field, data_found)
                call field_expression%set(permittivity, gap, size_x, surface_density)
            class default
                call error_exit("field_expression type unknown.")
        end select
    end subroutine set_field_expression

    subroutine destroy(field_expression)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression

        if (allocated(field_expression)) deallocate(field_expression)
    end subroutine destroy

end module procedures_field_expression_factory
