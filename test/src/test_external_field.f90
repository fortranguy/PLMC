module procedures_external_field_write

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_geometry, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use class_external_field, only: External_Field_Facade

implicit none

private
public write_field

contains

    subroutine write_field(field_unit, periodic_box, external_field, delta)
        integer, intent(in) :: field_unit
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(External_Field_Facade), intent(in) :: external_field
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
use module_data, only: test_file_exists, test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use types_field_parameters, only: Abstract_Field_Parameters, Constant_Field_Parameters
use class_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
                                       Concrete_Parallelepiped_Domain
use class_external_field, only: External_Field_Facade
use procedures_external_field_write, only: write_field

implicit none

    type(External_Field_Facade) :: external_field
    class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Field_Expression), allocatable :: field_expression
    class(Abstract_Field_Parameters), allocatable :: field_parameters
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field, field_name
    logical :: found
    real(DP), allocatable :: field_vector(:), box_size(:), domain_origin(:), domain_size(:), &
                             delta(:)
    integer :: field_unit

    call json_initialize()
    data_filename = "external_field.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    data_field = "External Field.name"
    call input_data%get(data_field, field_name, found)
    call test_data_found(data_field, found)

    select case (field_name)
        case ("constant")
            allocate(Constant_Field_Parameters :: field_parameters)
        case default
            call error_exit("Field not yet implemented?")
    end select
    select type (field_parameters)
        type is (Constant_Field_Parameters)
            data_field = "External Field.vector"
            call input_data%get(data_field, field_vector, found)
            call test_data_found(data_field, found)
            field_parameters%vector = field_vector
            deallocate(field_vector)
            allocate(Constant_Field_Expression :: field_expression)
        class default
            call error_exit("field_parameters unknown.")
    end select
    call field_expression%set(field_parameters)

    data_field = "Periodic Box.size"
    call input_data%get(data_field, box_size, found)
    call test_data_found(data_field, found)
    allocate(XYZ_Periodic_Box :: periodic_box)
    call periodic_box%set_size(box_size)
    deallocate(box_size)

    data_field = "Parallelepiped Domain.origin"
    call input_data%get(data_field, domain_origin, found)
    call test_data_found(data_field, found)
    data_field = "Parallelepiped Domain.size"
    call input_data%get(data_field, domain_size, found)
    call test_data_found(data_field, found)
    allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
    call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)

    data_field = "External Field.delta"
    call input_data%get(data_field, delta, found)
    call test_data_found(data_field, found)
    call external_field%construct(parallelepiped_domain, field_expression)
    open(newunit=field_unit, recl=4096, file="constant_field.out", action="write")
    call write_field(field_unit, periodic_box, external_field, delta)
    close(field_unit)

    call external_field%destroy()
    call parallelepiped_domain%destroy()
    deallocate(parallelepiped_domain)
    deallocate(periodic_box)
    deallocate(field_expression)
    deallocate(field_parameters)
    deallocate(field_name)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_external_field
