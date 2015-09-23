module procedures_box_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_errors, only: error_exit
use procedures_checks, only: check_3d_array
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_temperature, only: Abstract_Temperature, Concrete_Temperature
use types_field_parameters, only: Abstract_Field_Parameters, Constant_Field_Parameters, &
    Null_Field_Parameters
use class_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression, &
    Null_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Null_Parallelepiped_Domain, Box_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field, Concrete_External_Field, &
    Null_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice, &
    Null_Reciprocal_Lattice
use types_box, only: Box_Wrapper

implicit none

private
public :: box_factory_construct, box_factory_destroy

contains

    subroutine box_factory_construct(box, input_data, prefix)
        type(Box_Wrapper), intent(out) :: box
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_and_set_periodic_box(box%periodic_box, input_data, prefix)
        call allocate_temperature(box%temperature, input_data, prefix)
        call allocate_field_parameters_and_expression(box%field_parameters, box%field_expression, &
            input_data, prefix)
        call set_field_parameters_and_expression(box%field_expression, box%field_parameters, &
            input_data, prefix)
        call allocate_and_construct_parallelepiped_domain(box%parallelepiped_domain, &
            input_data, prefix//".External Field", box%periodic_box)
        call allocate_external_field(box%external_field, input_data, prefix)
        call box%external_field%construct(box%parallelepiped_domain, box%field_expression)
        call allocate_and_construct_reciprocal_lattice(box%reciprocal_lattice, input_data, prefix, &
            box%periodic_box)
    end subroutine box_factory_construct

    subroutine allocate_and_set_periodic_box(periodic_box, input_data, prefix)
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: box_periodicity
        real(DP), allocatable :: box_size(:)

        data_field = prefix//".Periodic Box.periodicity"
        call input_data%get(data_field, box_periodicity, data_found)
        call test_data_found(data_field, data_found)
        select case (box_periodicity)
            case ("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_box)
            case ("XY")
                allocate(XY_Periodic_Box :: periodic_box)
            case default
                call error_exit(data_field//" unknown. Choose between: 'XYZ' and 'XY'")
        end select
        deallocate(box_periodicity)
        data_field = prefix//".Periodic Box.size"
        call input_data%get(data_field, box_size, data_found)
        call test_data_found(data_field, data_found)
        deallocate(data_field)
        call periodic_box%set(box_size)
    end subroutine allocate_and_set_periodic_box

    subroutine allocate_temperature(temperature, input_data, prefix)
        class(Abstract_Temperature), allocatable, intent(out) :: temperature
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: temperature_value

        data_field = prefix//".Temperature"
        call input_data%get(data_field, temperature_value, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
        deallocate(data_field)
    end subroutine allocate_temperature

    subroutine allocate_field_parameters_and_expression(field_parameters, field_expression, &
        input_data, prefix)
        class(Abstract_Field_Parameters), allocatable, intent(out) :: field_parameters
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, field_name
        logical :: data_found

        if (apply_external_field(input_data, prefix)) then
            data_field = prefix//".External Field.name"
            call input_data%get(data_field, field_name, data_found)
            call test_data_found(data_field, data_found)
            select case (field_name)
                case ("constant")
                    allocate(Constant_Field_Parameters :: field_parameters)
                    allocate(Constant_Field_Expression :: field_expression)
                case default
                    call error_exit(field_name//" field_name unknown. Choose: 'constant'.")
            end select
            deallocate(field_name)
            deallocate(data_field)
        else
            allocate(Null_Field_Parameters :: field_parameters)
            allocate(Null_Field_Expression :: field_expression)
        end if
    end subroutine allocate_field_parameters_and_expression

    subroutine set_field_parameters_and_expression(field_expression, field_parameters, input_data, &
        prefix)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression
        class(Abstract_Field_Parameters), intent(inout) :: field_parameters
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: field_vector(:)

        if (apply_external_field(input_data, prefix)) then
            select type (field_parameters)
                type is (Constant_Field_Parameters)
                    data_field = prefix//".External Field.vector"
                    call input_data%get(data_field, field_vector, data_found)
                    call test_data_found(data_field, data_found)
                    call check_3d_array("Concrete_Box_Factory", "field_vector", field_vector)
                    field_parameters%vector = field_vector
                    deallocate(field_vector)
                    deallocate(data_field)
            end select
        end if
        call field_expression%set(field_parameters)
    end subroutine set_field_parameters_and_expression

    subroutine allocate_and_construct_parallelepiped_domain(parallelepiped_domain, input_data, &
        prefix, periodic_box)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: domain_name
        real(DP), allocatable :: domain_origin(:), domain_size(:)

        data_field = prefix//".Parallelepiped Domain.name"
        call input_data%get(data_field, domain_name, data_found)
        call test_data_found(data_field, data_found)
        select case(domain_name)
            case ("domain")
                data_field = prefix//".Parallelepiped Domain.origin"
                call input_data%get(data_field, domain_origin, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//".Parallelepiped Domain.size"
                call input_data%get(data_field, domain_size, data_found)
                call test_data_found(data_field, data_found)
                allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
            case ("box")
                allocate(Box_Parallelepiped_Domain :: parallelepiped_domain)
            case ("null")
                allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
            case default
                call error_exit(domain_name//" domain_name unknown."//&
                    "Choose among: 'domain', 'box', 'null'.")
        end select
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
        if (allocated(domain_size)) deallocate(domain_size)
        if (allocated(domain_origin)) deallocate(domain_origin)
        deallocate(domain_name)
        deallocate(data_field)
    end subroutine allocate_and_construct_parallelepiped_domain

    subroutine allocate_external_field(external_field, input_data, prefix)
        class(Abstract_External_Field), allocatable, intent(out) :: external_field
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (apply_external_field(input_data, prefix)) then
            allocate(Concrete_External_Field :: external_field)
        else
            allocate(Null_External_Field :: external_field)
        end if
    end subroutine allocate_external_field

    function apply_external_field(input_data, prefix)
        logical :: apply_external_field
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//".External Field.apply"
        call input_data%get(data_field, apply_external_field, data_found)
        call test_data_found(data_field, data_found)
        deallocate(data_field)
    end function apply_external_field

    subroutine allocate_and_construct_reciprocal_lattice(reciprocal_lattice, input_data, prefix, &
        periodic_box)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(out) :: reciprocal_lattice
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer, allocatable :: numbers(:)

        if (use_reciprocal_lattice(input_data, prefix)) then
            data_field = prefix//".Reciprocal Lattice.numbers"
            call input_data%get(data_field, numbers, data_found)
            call test_data_found(data_field, data_found)
            deallocate(data_field)
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
        end if
        call reciprocal_lattice%construct(periodic_box, numbers)
        if (allocated(numbers)) deallocate(numbers)
    end subroutine allocate_and_construct_reciprocal_lattice

    function use_reciprocal_lattice(input_data, prefix)
        logical :: use_reciprocal_lattice
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//".Reciprocal Lattice.use"
        call input_data%get(data_field, use_reciprocal_lattice, data_found)
        call test_data_found(data_field, data_found)
        deallocate(data_field)
    end function use_reciprocal_lattice

    subroutine box_factory_destroy(box)
        type(Box_Wrapper), intent(inout) :: box

        call box%reciprocal_lattice%destroy()
        if (allocated(box%reciprocal_lattice)) deallocate(box%reciprocal_lattice)
        call box%external_field%destroy()
        if (allocated(box%external_field)) deallocate(box%external_field)
        call box%parallelepiped_domain%destroy()
        if (allocated(box%parallelepiped_domain)) deallocate(box%parallelepiped_domain)
        if (allocated(box%field_expression)) deallocate(box%field_expression)
        if (allocated(box%field_parameters)) deallocate(box%field_parameters)
        if (allocated(box%temperature)) deallocate(box%temperature)
        if (allocated(box%periodic_box)) deallocate(box%periodic_box)
    end subroutine box_factory_destroy

end module procedures_box_factory
