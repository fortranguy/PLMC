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
use procedures_checks, only: check_file_exists, check_data_found
use procedures_errors, only: error_exit
use class_particles_diameter, only: Abstract_Particles_Diameter
use procedures_particles_factory, only: particles_factory_create, particles_factory_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use procedures_potential_expression_write, only: write_potential_expression

implicit none

    class(Abstract_Particles_Diameter), allocatable :: diameter
    class(Abstract_Potential_Expression), allocatable :: potential_expression
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    real(DP) :: min_distance, max_distance, delta_distance
    integer :: potential_unit

    call json_initialize()
    data_filename = "potential_expression.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call particles_factory_create(diameter, .true., input_data, "Test Potential Expression")
    call short_potential_factory_create(potential_expression, .true., input_data, &
        "Test Potential Expression")

    data_field = "Test Potential Expression.Potential.minimum distance"
    call input_data%get(data_field, min_distance, data_found)
    call check_data_found(data_field, data_found)
    data_field = "Test Potential Expression.Potential.maximum distance"
    call input_data%get(data_field, max_distance, data_found)
    call check_data_found(data_field, data_found)
    data_field = "Test Potential Expression.Potential.delta distance"
    call input_data%get(data_field, delta_distance, data_found)
    call check_data_found(data_field, data_found)

    open(newunit=potential_unit, recl=4096, file="potential_expression.out", action="write")
    call write_potential_expression(potential_expression, min_distance, max_distance, &
                                    delta_distance, potential_unit)
    close(potential_unit)

    call particles_factory_destroy(diameter)
    call short_potential_factory_destroy(potential_expression)
    deallocate(data_field)
    call input_data%destroy()

end program test_potential_expression
