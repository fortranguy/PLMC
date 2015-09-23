module class_box_factory

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
    Concrete_Parallelepiped_Domain, Null_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field, Concrete_External_Field, &
    Null_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice, &
    Null_Reciprocal_Lattice
use types_box, only: Box_Wrapper

implicit none

private

    type, public :: Concrete_Box_Factory
    private
        type(json_file), pointer :: input_data
        character(len=:), allocatable :: prefix
        logical :: use_external_field
        logical :: use_reciprocal_lattice
    contains
        procedure :: allocate => Concrete_Box_Factory_allocate
        procedure, private :: allocate_periodic_box => Concrete_Box_Factory_allocate_periodic_box
        procedure, private :: allocate_temperature => Concrete_Box_Factory_allocate_temperature
        procedure, private :: allocate_field_parameters => &
            Concrete_Box_Factory_allocate_field_parameters
        procedure, private :: allocate_field_expression => &
            Concrete_Box_Factory_allocate_field_expression
        procedure, private :: allocate_parallelepiped_domain => &
            Concrete_Box_Factory_allocate_parallelepiped_domain
        procedure, private :: allocate_external_field => &
            Concrete_Box_Factory_allocate_external_field
        procedure, private :: allocate_reciprocal_lattice => &
            Concrete_Box_Factory_allocate_reciprocal_lattice
        procedure :: construct => Concrete_Box_Factory_construct
        procedure, private :: construct_parallelepiped_domain => &
            Concrete_Box_Factory_construct_parallelepiped_domain
        procedure, private :: construct_reciprocal_lattice => &
            Concrete_Box_Factory_construct_reciprocal_lattice
        procedure :: destroy => Concrete_Box_Factory_destroy
    end type Concrete_Box_Factory

contains

    subroutine Concrete_Box_Factory_allocate(this, box, input_data, prefix)
        class(Concrete_Box_Factory), intent(out) :: this
        type(Box_Wrapper), intent(out) :: box
        type(json_file), target, intent(in) :: input_data
        character(len=*), intent(in) :: prefix

        this%input_data => input_data
        this%prefix = prefix
        call this%allocate_periodic_box(box%periodic_box)
        call this%allocate_temperature(box%temperature)
        call this%allocate_field_parameters(box%field_parameters)
        call this%allocate_field_expression(box%field_expression, box%field_parameters)
        call this%allocate_parallelepiped_domain(box%parallelepiped_domain)
        call this%allocate_external_field(box%external_field)
        call this%allocate_reciprocal_lattice(box%reciprocal_lattice)
    end subroutine Concrete_Box_Factory_allocate

    subroutine Concrete_Box_Factory_allocate_periodic_box(this, periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: box_name
        real(DP), allocatable :: box_size(:)

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
        data_field = this%prefix//".Periodic Box.size"
        call this%input_data%get(data_field, box_size, data_found)
        call test_data_found(data_field, data_found)
        call periodic_box%set(box_size)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_periodic_box

    subroutine Concrete_Box_Factory_allocate_temperature(this, temperature)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Temperature), allocatable, intent(out) :: temperature

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: temperature_value

        data_field = this%prefix//".Temperature"
        call this%input_data%get(data_field, temperature_value, data_found)
        call test_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_temperature

    subroutine Concrete_Box_Factory_allocate_field_parameters(this, field_parameters)
        class(Concrete_Box_Factory), intent(inout) :: this
        class(Abstract_Field_Parameters), allocatable, intent(out) :: field_parameters

        character(len=:), allocatable :: data_field, field_name
        logical :: data_found

        data_field = this%prefix//".Field.name"
        call this%input_data%get(data_field, field_name, data_found)
        call test_data_found(data_field, data_found)
        this%use_external_field = .true.
        select case (field_name)
            case("constant")
                allocate(Constant_Field_Parameters :: field_parameters)
            case("null")
                this%use_external_field = .false.
                allocate(Null_Field_Parameters :: field_parameters)
            case default
                call error_exit(field_name//"unknown. Choose between: 'constant' and 'null'.")
        end select
        deallocate(field_name)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_field_parameters

    subroutine Concrete_Box_Factory_allocate_field_expression(this, field_expression, &
        field_parameters)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        class(Abstract_Field_Parameters), intent(inout) :: field_parameters

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: field_vector(:)

        if (this%use_external_field) then
            select type (field_parameters)
                type is (Constant_Field_Parameters)
                    data_field = this%prefix//".Field.vector"
                    call this%input_data%get(data_field, field_vector, data_found)
                    call test_data_found(data_field, data_found)
                    call check_3d_array("Concrete_Box_Factory", "field_vector", field_vector)
                    field_parameters%vector = field_vector
                    deallocate(field_vector)
                    allocate(Constant_Field_Expression :: field_expression)
            end select
        else
            allocate(Null_Field_Expression :: field_expression)
        end if
        call field_expression%set(field_parameters)
        if (allocated(data_field)) deallocate(data_field)
    end subroutine Concrete_Box_Factory_allocate_field_expression

    subroutine Concrete_Box_Factory_allocate_parallelepiped_domain(this, parallelepiped_domain)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain

        if (this%use_external_field) then
            allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if
    end subroutine Concrete_Box_Factory_allocate_parallelepiped_domain

    subroutine Concrete_Box_Factory_allocate_external_field(this, external_field)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_External_Field), allocatable, intent(out) :: external_field

        if (this%use_external_field) then
            allocate(Concrete_External_Field :: external_field)
        else
            allocate(Null_External_Field :: external_field)
        end if
    end subroutine Concrete_Box_Factory_allocate_external_field

    subroutine Concrete_Box_Factory_allocate_reciprocal_lattice(this, reciprocal_lattice)
        class(Concrete_Box_Factory), intent(inout) :: this
        class(Abstract_Reciprocal_Lattice), allocatable, intent(out) :: reciprocal_lattice

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = this%prefix//".Reciprocal Lattice.use"
        call this%input_data%get(data_field, this%use_reciprocal_lattice, data_found)
        call test_data_found(data_field, data_found)
        if (this%use_reciprocal_lattice) then
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
        end if
    end subroutine Concrete_Box_Factory_allocate_reciprocal_lattice

    subroutine Concrete_Box_Factory_construct(this, box)
        class(Concrete_Box_Factory), intent(in) :: this
        type(Box_Wrapper), intent(inout) :: box

        call this%construct_parallelepiped_domain(box%parallelepiped_domain, box%periodic_box)
        call box%external_field%construct(box%parallelepiped_domain, box%field_expression)
        call this%construct_reciprocal_lattice(box%reciprocal_lattice, box%periodic_box)
    end subroutine Concrete_Box_Factory_construct

    subroutine Concrete_Box_Factory_construct_parallelepiped_domain(this, parallelepiped_domain, &
        periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Parallelepiped_Domain), intent(inout) :: parallelepiped_domain
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: domain_origin(:), domain_size(:)

        if (.not. this%use_external_field) return
        data_field = this%prefix//".Field.Parallelepiped Domain.origin"
        call this%input_data%get(data_field, domain_origin, data_found)
        call test_data_found(data_field, data_found)
        data_field = this%prefix//".Field.Parallelepiped Domain.size"
        call this%input_data%get(data_field, domain_size, data_found)
        call test_data_found(data_field, data_found)
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
        deallocate(domain_size)
        deallocate(domain_origin)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_construct_parallelepiped_domain

    subroutine Concrete_Box_Factory_construct_reciprocal_lattice(this, reciprocal_lattice, &
        periodic_box)
        class(Concrete_Box_Factory), intent(in) :: this
        class(Abstract_Reciprocal_Lattice), intent(inout) :: reciprocal_lattice
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer, allocatable :: numbers(:)

        if (.not. this%use_reciprocal_lattice) return
        data_field = this%prefix//".Reciprocal Lattice.numbers"
        call this%input_data%get(data_field, numbers, data_found)
        call test_data_found(data_field, data_found)
        call reciprocal_lattice%construct(periodic_box, numbers)
        deallocate(numbers)
        deallocate(data_field)
    end subroutine Concrete_Box_Factory_construct_reciprocal_lattice

    subroutine Concrete_Box_Factory_destroy(this, box)
        class(Concrete_Box_Factory), intent(inout) :: this
        type(Box_Wrapper), intent(inout) :: box

        if (allocated(box%reciprocal_lattice)) deallocate(box%reciprocal_lattice)
        if (allocated(box%external_field)) deallocate(box%external_field)
        if (allocated(box%parallelepiped_domain)) deallocate(box%parallelepiped_domain)
        if (allocated(box%field_expression)) deallocate(box%field_expression)
        if (allocated(box%field_parameters)) deallocate(box%field_parameters)
        if (allocated(box%temperature)) deallocate(box%temperature)
        if (allocated(box%periodic_box)) deallocate(box%periodic_box)
        if (allocated(this%prefix)) deallocate(this%prefix)
        this%input_data => null()
    end subroutine Concrete_Box_Factory_destroy

end module class_box_factory
