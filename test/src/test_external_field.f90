module procedures_external_field_write

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use class_external_field, only: Abstract_External_Field

implicit none

private
public write_field

contains

    subroutine write_field(field_unit, periodic_box, external_field, delta)
        integer, intent(in) :: field_unit
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_External_Field), intent(in) :: external_field
        real(DP), intent(in) :: delta(:)

        integer :: indices(num_dimensions), i, j, k
        real(DP) :: position(num_dimensions)

        indices = int(periodic_box%get_size()/delta)/2
        do k=-indices(3), indices(3)
        do j=-indices(2), indices(2)
        do i=-indices(1), indices(1)
            position = [i, j, k] * delta
            write(field_unit, *) position, external_field%get(position)
        end do
        end do
        end do
    end subroutine write_field

end module procedures_external_field_write

program test_external_field

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use procedures_checks, only: check_file_exists, check_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_field_expression, only: Abstract_Field_Expression
use procedures_environment_factory, only: environment_factory_create, environment_factory_destroy
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
                                       Concrete_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field, Concrete_External_Field
use procedures_external_field_write, only: write_field

implicit none

    class(Abstract_External_Field), allocatable :: external_field
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Field_Expression), allocatable :: field_expression
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    real(DP), allocatable :: delta(:)
    integer :: field_unit

    call json_initialize()
    data_filename = "external_field.json"
    call check_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    call environment_factory_create(periodic_box, input_data, "Environment.")
    call environment_factory_create(field_expression, input_data, "Environment.")
    call environment_factory_create(parallelepiped_domain, input_data, &
        "Environment.External Field.", periodic_box)

    call environment_factory_create(external_field, input_data, "Environment.", &
        parallelepiped_domain, field_expression)
    call environment_factory_destroy(field_expression)
    open(newunit=field_unit, recl=4096, file="external_field.out", action="write")
    data_field = "Environment.External Field.delta"
    call input_data%get(data_field, delta, data_found)
    call check_data_found(data_field, data_found)
    call write_field(field_unit, periodic_box, external_field, delta)
    deallocate(delta)
    close(field_unit)

    call environment_factory_destroy(external_field)
    call environment_factory_destroy(parallelepiped_domain)
    call environment_factory_destroy(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_external_field
