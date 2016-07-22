module procedures_potential_expression_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_potential_expression, only: Abstract_Potential_Expression, Lennard_Jones_Expression, &
    Null_Potential_Expression

implicit none

private
public :: create, destroy

contains


    subroutine create(expression, interact, generating_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        call allocate(expression, interact, generating_data, prefix)
        call set(expression, generating_data, prefix)
    end subroutine create

    subroutine allocate(expression, interact, generating_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        if (interact) then
            data_field = prefix//"name"
            call generating_data%get(data_field, potential_name, data_found)
            call check_data_found(data_field, data_found)
            select case (potential_name)
                case ("null")
                    allocate(Null_Potential_Expression :: expression)
                case ("LJ")
                    allocate(Lennard_Jones_Expression :: expression)
                case default
                    call error_exit(potential_name//" unknown potential_name. "//&
                        "Choose between: 'null' and LJ.")
            end select
        else
            allocate(Null_Potential_Expression :: expression)
        end if
    end subroutine allocate

    subroutine set(expression, generating_data, prefix)
        class(Abstract_Potential_Expression), intent(inout) :: expression
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        select type(expression)
            type is (Lennard_Jones_Expression)
                block
                    real(DP) :: epsilon, sigma
                    data_field = prefix//"epsilon"
                    call generating_data%get(data_field, epsilon, data_found)
                    call check_data_found(data_field, data_found)
                    data_field = prefix//"sigma"
                    call generating_data%get(data_field, sigma, data_found)
                    call check_data_found(data_field, data_found)
                    call expression%set(epsilon, sigma)
                end block
            type is (Null_Potential_Expression)
            class default
                call error_exit("procedures_potential_expression_factory: expression type unknown.")
        end select
    end subroutine set

    subroutine destroy(expression)
        class(Abstract_Potential_Expression), allocatable, intent(inout) :: expression

        if (allocated(expression)) deallocate(expression)
    end subroutine destroy

end module procedures_potential_expression_factory
