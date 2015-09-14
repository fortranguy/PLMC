module procedures_potential_expression_write

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_potential_expression, only: Abstract_Potential_Expression

implicit none

private
public write_potential_expression

contains

    subroutine write_potential_expression(potential_expression, min_distance, max_distance, &
                                          delta_distance, write_unit)
        class(Abstract_Potential_Expression), intent(in) :: potential_expression
        real(DP), intent(in) :: min_distance, max_distance, delta_distance
        integer, intent(in) :: write_unit

        real(DP) :: distance

        distance = min_distance
        do while(distance < max_distance)
            write(write_unit, *) distance, potential_expression%get(distance)
            distance = distance + delta_distance
        end do
    end subroutine write_potential_expression

end module procedures_potential_expression_write

program test_potential_expression

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use types_potential_parameters, only: Abstract_Potential_Parameters, Null_Potential_Paramters, &
    Lennard_Jones_Parameters
use class_potential_expression, only: Abstract_Potential_Expression, Null_Potential_Expression, &
    Lennard_Jones_Expression
use procedures_potential_expression_write, only: write_potential_expression

implicit none

    class(Abstract_Potential_Expression), allocatable :: potential_expression
    class(Abstract_Potential_Parameters), allocatable :: potential_paramters
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field, potential_name
    logical :: data_found

    real(DP) :: min_distance, max_distance, delta_distance
    integer :: potential_unit

    call json_initialize()
    data_filename = "potential_expression.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    data_field = "Potential.name"
    call input_data%get(data_field, potential_name, data_found)
    call test_data_found(data_field, data_found)

    select case(potential_name)
        case ("null")
            allocate(Null_Potential_Paramters :: potential_paramters)
        case ("LJ")
            allocate(Lennard_Jones_Parameters :: potential_paramters)
        case default
            call error_exit(data_field//" unkown.")
    end select

    select type(potential_paramters)
        type is (Null_Potential_Paramters)
            allocate(Null_Potential_Expression :: potential_expression)
        type is (Lennard_Jones_Parameters)
            data_field = "Potential.epsilon"
            call input_data%get(data_field, potential_paramters%epsilon, data_found)
            call test_data_found(data_field, data_found)
            data_field = "Potential.sigma"
            call input_data%get(data_field, potential_paramters%sigma, data_found)
            call test_data_found(data_field, data_found)
            allocate(Lennard_Jones_Expression :: potential_expression)
    end select
    call potential_expression%set(potential_paramters)

    data_field = "Potential.minimum distance"
    call input_data%get(data_field, min_distance, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Potential.maximum distance"
    call input_data%get(data_field, max_distance, data_found)
    call test_data_found(data_field, data_found)
    data_field = "Potential.delta distance"
    call input_data%get(data_field, delta_distance, data_found)
    call test_data_found(data_field, data_found)

    open(newunit=potential_unit, recl=4096, file="potential_expression.out", action="write")
    call write_potential_expression(potential_expression, min_distance, max_distance, &
                                    delta_distance, potential_unit)
    close(potential_unit)

    deallocate(potential_name)
    deallocate(potential_expression)
    deallocate(potential_paramters)
    deallocate(data_field)
    call input_data%destroy()

end program test_potential_expression
