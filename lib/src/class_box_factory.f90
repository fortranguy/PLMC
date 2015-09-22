module class_box_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_temperature, only: Abstract_Temperature, Concrete_Temperature
use types_field_parameters, only: Abstract_Field_Parameters, Constant_Field_Parameters, &
    Null_Field_Parameters
use class_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression, &
    Null_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain
use types_box, only: Box_Wrapper

implicit none

private

    type, public :: Concrete_Box_Factory
    private
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
    contains
        procedure :: allocate => Concrete_Box_Factory_allocate
        procedure :: allocate_periodic_box => Concrete_Box_Factory_allocate_periodic_box
        procedure :: allocate_temperature => Concrete_Box_Factory_allocate_temperature
        procedure :: allocate_field_expression => Concrete_Box_Factory_allocate_field_expression
        !procedure :: construct => Concrete_Box_Factory_construct
        procedure :: destroy => Concrete_Box_Factory_destroy
    end type Concrete_Box_Factory

contains

    subroutine Concrete_Box_Factory_allocate(this, box, input_data, prefix)
        class(Concrete_Box_Factory), intent(out) :: this
        type(Box_Wrapper), intent(out) :: box
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        this%input_data = input_data
        this%prefix = prefix
        call this%allocate_periodic_box(box%periodic_box)
    end subroutine Concrete_Box_Factory_allocate

    subroutine Concrete_Box_Factory_allocate_periodic_box(this, periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: box_name

        data_field = this%prefix//".Periodic Box.name"
        call this%input_data%get(data_field, box_name, data_found)
        call test_data_found(data_field, data_found)
        select case(box_name)
            case("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_box)
            case("XY")
                allocate(XY_Periodic_Box :: periodic_box)
            case default
                call error_exit(data_field//" unkown. Choose between: 'XYZ' and 'XY'")
        end select
        deallocate(box_name)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_periodic_box

    subroutine Concrete_Box_Factory_allocate_field_expression(this, field_expression)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression

        class(Abstract_Field_Parameters), allocatable :: field_parameters
        character(len=:), allocatable :: data_field, field_name
        logical :: data_found
        real(DP), allocatable :: field_vector(:)

        data_field = this%prefix//".Field.name"
        call this%input_data%get(data_field, field_name, data_found)
        call test_data_found(data_field, data_found)

        select case (field_name)
            case("constant")
                allocate(Constant_Field_Parameters :: field_parameters)
            case("null")
                allocate(Null_Field_Parameters :: field_parameters)
            case default
                call error_exit(field_name//"unknown. Choose between: 'constant' and 'null'.")
        end select
        select type (field_parameters)
            type is (Constant_Field_Parameters)
                data_field = this%prefix//".Field.vector"
                call this%input_data%get(data_field, field_vector, data_found)
                call test_data_found(data_field, data_found)
                field_parameters%vector = field_vector
                deallocate(field_vector)
                allocate(Constant_Field_Expression :: field_expression)
            type is (Null_Field_Parameters)
                allocate(Null_Field_Expression :: field_expression)
            class default
                call error_exit("Field parameters unknown.")
        end select
        call field_expression%set(field_parameters)
        deallocate(field_parameters)
        deallocate(field_name)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_field_expression

    subroutine Concrete_Box_Factory_allocate_parallelepiped_domain(this, parallelepiped_domain)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain

        allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
    end subroutine Concrete_Box_Factory_allocate_parallelepiped_domain

    subroutine Concrete_Box_Factory_allocate_temperature(this, temperature)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Temperature), allocatable, intent(out) :: temperature

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: temperature_value

        data_field = this%prefix//"Temperature"
        call this%input_data%get(data_field, temperature_value, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_temperature

    subroutine Concrete_Box_Factory_construct_parallelepiped_domain(this, parallelepiped_domain, &
        periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Parallelepiped_Domain), intent(inout) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: domain_origin(:), domain_size(:)

        data_field = this%prefix//".Parallelepiped Domain.origin"
        call this%input_data%get(data_field, domain_origin, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Parallelepiped Domain.size"
        call this%input_data%get(data_field, domain_size, data_found)
        call test_data_found(data_field, data_found)
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
        deallocate(domain_origin)
        deallocate(domain_size)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_construct_parallelepiped_domain

    subroutine Concrete_Box_Factory_destroy(this, box)
        class(Concrete_Box_Factory), intent(inout) :: this
        type(Box_Wrapper), intent(inout) :: box

        if (allocated(box%temperature)) deallocate(box%temperature)
        if (allocated(box%periodic_box)) deallocate(box%periodic_box)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Concrete_Box_Factory_destroy

end module class_box_factory
