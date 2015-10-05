module procedures_environment_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found, check_3d_array
use class_periodic_box, only: Abstract_Periodic_Box, &
    XYZ_Periodic_Box, XY_Periodic_Box
use class_temperature, only: Abstract_Temperature, &
    Concrete_Temperature
use class_field_expression, only: Abstract_Field_Expression, &
    Constant_Field_Expression, Null_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Null_Parallelepiped_Domain, Box_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field, &
    Concrete_External_Field, Null_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, &
    Concrete_Reciprocal_Lattice, Null_Reciprocal_Lattice
use class_floor_penetration, only: Abstract_Floor_Penetration, &
    Flat_Floor_Penetration, Null_Floor_Penetration
use procedures_particles_factory, only: particles_factory_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use procedures_short_potential_factory, only: short_potential_factory_create, &
    short_potential_factory_destroy
use class_walls_potential, only: Abstract_Walls_Potential, &
    Concrete_Walls_Potential, Null_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use procedures_property_inquirers, only: apply_external_field, use_reciprocal_lattice, use_walls

implicit none

private
public :: environment_factory_create, environment_factory_destroy

interface environment_factory_create
    module procedure :: environment_factory_create_all
    module procedure :: allocate_and_set_periodic_box
    module procedure :: allocate_and_set_temperature
    module procedure :: allocate_and_set_field_expression
    module procedure :: allocate_and_construct_parallelepiped_domain
    module procedure :: allocate_and_construct_external_field
    module procedure :: allocate_and_construct_reciprocal_lattice
    module procedure :: allocate_and_set_floor_penetration
    module procedure :: allocate_and_construct_walls_potential
end interface environment_factory_create

interface environment_factory_destroy
    module procedure :: destroy_and_deallocate_walls_potential
    module procedure :: deallocate_floor_penetration
    module procedure :: destroy_and_deallocate_reciprocal_lattice
    module procedure :: destroy_and_deallocate_external_field
    module procedure :: destroy_and_deallocate_parallelepiped_domain
    module procedure :: deallocate_field_expression
    module procedure :: deallocate_temperature
    module procedure :: deallocate_periodic_box
    module procedure :: environment_factory_destroy_all
end interface environment_factory_destroy

contains

    subroutine environment_factory_create_all(environment, input_data, prefix)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        logical :: apply_field

        call environment_factory_create(environment%periodic_box, input_data, prefix)
        call environment_factory_create(environment%temperature, input_data, prefix)
        apply_field = apply_external_field(input_data, prefix)
        call environment_factory_create(field_expression, apply_field, input_data, prefix)
        call environment_factory_create(parallelepiped_domain, apply_field, &
            environment%periodic_box, input_data, prefix//"External Field.")
        call environment_factory_create(environment%external_field, apply_field, &
            parallelepiped_domain, field_expression)
        call environment_factory_destroy(field_expression)
        call environment_factory_destroy(parallelepiped_domain)
        call environment_factory_create(environment%reciprocal_lattice, input_data, prefix, &
            environment%periodic_box)
        call environment_factory_create(floor_penetration, input_data, prefix)
        call environment_factory_create(environment%walls_potential, input_data, prefix, &
            environment%periodic_box, floor_penetration)
        call environment_factory_destroy(floor_penetration)
    end subroutine environment_factory_create_all

    subroutine allocate_and_set_periodic_box(periodic_box, input_data, prefix)
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: box_periodicity
        real(DP), allocatable :: box_size(:)

        data_field = prefix//"Box.periodicity"
        call input_data%get(data_field, box_periodicity, data_found)
        call check_data_found(data_field, data_found)
        select case (box_periodicity)
            case ("XYZ")
                allocate(XYZ_Periodic_Box :: periodic_box)
            case ("XY")
                allocate(XY_Periodic_Box :: periodic_box)
            case default
                call error_exit(data_field//" unknown. Choose between: 'XYZ' and 'XY'")
        end select
        deallocate(box_periodicity)
        data_field = prefix//"Box.size"
        call input_data%get(data_field, box_size, data_found)
        call check_data_found(data_field, data_found)
        deallocate(data_field)
        call periodic_box%set(box_size)
    end subroutine allocate_and_set_periodic_box

    subroutine allocate_and_set_temperature(temperature, input_data, prefix)
        class(Abstract_Temperature), allocatable, intent(out) :: temperature
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: temperature_value

        data_field = prefix//"Thermostat.temperature"
        call input_data%get(data_field, temperature_value, data_found)
        call check_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
        deallocate(data_field)
    end subroutine allocate_and_set_temperature

    subroutine allocate_and_set_field_expression(field_expression, apply_field, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        logical, intent(in) :: apply_field
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_field_expression(field_expression, apply_field, input_data, prefix)
        call set_field_expression(field_expression, input_data, prefix)
    end subroutine allocate_and_set_field_expression

    subroutine allocate_field_expression(field_expression, apply_field, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        logical, intent(in) :: apply_field
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, field_name
        logical :: data_found
        if (apply_field) then
            data_field = prefix//"External Field.name"
            call input_data%get(data_field, field_name, data_found)
            call check_data_found(data_field, data_found)
            select case (field_name)
                case ("constant")
                    allocate(Constant_Field_Expression :: field_expression)
                case default
                    call error_exit(field_name//" field_name unknown. Choose: 'constant'.")
            end select
            deallocate(field_name)
            deallocate(data_field)
        else
            allocate(Null_Field_Expression :: field_expression)
        end if
    end subroutine allocate_field_expression

    subroutine set_field_expression(field_expression, input_data, prefix)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP), allocatable :: field_vector(:)

        select type (field_expression)
            type is (Null_Field_Expression)
                call field_expression%set()
            type is (Constant_Field_Expression)
                data_field = prefix//"External Field.vector"
                call input_data%get(data_field, field_vector, data_found)
                call check_data_found(data_field, data_found)
                call field_expression%set(field_vector)
                deallocate(field_vector)
        end select
        if (allocated(data_field)) deallocate(data_field)
    end subroutine set_field_expression

    subroutine allocate_and_construct_parallelepiped_domain(parallelepiped_domain, needed, &
        periodic_box, input_data, prefix)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        logical, intent(in) :: needed
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: domain_name
        real(DP), allocatable :: domain_origin(:), domain_size(:)

        if (needed) then
            data_field = prefix//"Parallelepiped Domain.name"
            call input_data%get(data_field, domain_name, data_found)
            call check_data_found(data_field, data_found)
            select case(domain_name)
                case ("domain")
                    data_field = prefix//"Parallelepiped Domain.origin"
                    call input_data%get(data_field, domain_origin, data_found)
                    call check_data_found(data_field, data_found)
                    data_field = prefix//"Parallelepiped Domain.size"
                    call input_data%get(data_field, domain_size, data_found)
                    call check_data_found(data_field, data_found)
                    allocate(Concrete_Parallelepiped_Domain :: parallelepiped_domain)
                case ("box")
                    allocate(Box_Parallelepiped_Domain :: parallelepiped_domain)
                case default
                    call error_exit(domain_name//" domain_name unknown."//&
                        "Choose among: 'domain', 'box', 'null'.")
            end select
            deallocate(domain_name)
            deallocate(data_field)
        else
            domain_size = [0._DP, 0._DP, 0._DP]
            domain_origin = [0._DP, 0._DP, 0._DP]
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
        if (allocated(domain_size)) deallocate(domain_size)
        if (allocated(domain_origin)) deallocate(domain_origin)
    end subroutine allocate_and_construct_parallelepiped_domain

    subroutine allocate_and_construct_external_field(external_field, apply_field, &
        parallelepiped_domain, field_expression)
        class(Abstract_External_Field), allocatable, intent(out) :: external_field
        logical, intent(in) :: apply_field
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression

        if (apply_field) then
            allocate(Concrete_External_Field :: external_field)
        else
            allocate(Null_External_Field :: external_field)
        end if
        call external_field%construct(parallelepiped_domain, field_expression)
    end subroutine allocate_and_construct_external_field

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
            data_field = prefix//"Reciprocal Lattice.numbers"
            call input_data%get(data_field, numbers, data_found)
            call check_data_found(data_field, data_found)
            deallocate(data_field)
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
        end if
        call reciprocal_lattice%construct(periodic_box, numbers)
        if (allocated(numbers)) deallocate(numbers)
    end subroutine allocate_and_construct_reciprocal_lattice

    subroutine allocate_and_set_floor_penetration(floor_penetration, input_data, prefix)
        class(Abstract_Floor_Penetration), allocatable, intent(out) :: floor_penetration
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        character(len=:), allocatable :: walls_name

        if (use_walls(input_data, prefix)) then
            data_field = prefix//"Walls.name"
            call input_data%get(data_field, walls_name, data_found)
            call check_data_found(data_field, data_found)
            select case(walls_name)
                case ("flat")
                    allocate(Flat_Floor_Penetration :: floor_penetration)
                case default
                    call error_exit(walls_name//" walls_name unknown. Choose: 'flat'.")
            end select
        else
            allocate(Null_Floor_Penetration :: floor_penetration)
        end if
    end subroutine allocate_and_set_floor_penetration

    subroutine allocate_and_construct_walls_potential(walls_potential, input_data, prefix, &
        periodic_box, floor_penetration)
        class(Abstract_Walls_Potential), allocatable, intent(out) :: walls_potential
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: gap

        if (use_walls(input_data, prefix)) then
            allocate(Concrete_Walls_Potential :: walls_potential)
            data_field = prefix//"Walls.gap"
            call input_data%get(data_field, gap, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Walls_Potential :: walls_potential)
        end if
        call walls_potential%construct(periodic_box, gap, floor_penetration)
    end subroutine allocate_and_construct_walls_potential

    subroutine environment_factory_destroy_all(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call environment_factory_destroy(environment%walls_potential)
        call environment_factory_destroy(environment%reciprocal_lattice)
        call environment_factory_destroy(environment%external_field)
        call environment_factory_destroy(environment%temperature)
        call environment_factory_destroy(environment%periodic_box)
    end subroutine environment_factory_destroy_all

    subroutine deallocate_periodic_box(periodic_box)
        class(Abstract_Periodic_Box), allocatable, intent(inout) :: periodic_box

        if (allocated(periodic_box)) deallocate(periodic_box)
    end subroutine deallocate_periodic_box

    subroutine deallocate_temperature(temperature)
        class(Abstract_Temperature), allocatable, intent(inout) :: temperature

        if (allocated(temperature)) deallocate(temperature)
    end subroutine deallocate_temperature

    subroutine deallocate_field_expression(field_expression)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression

        if (allocated(field_expression)) deallocate(field_expression)
    end subroutine deallocate_field_expression

    subroutine destroy_and_deallocate_parallelepiped_domain(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: parallelepiped_domain

        call parallelepiped_domain%destroy()
        if (allocated(parallelepiped_domain)) deallocate(parallelepiped_domain)
    end subroutine destroy_and_deallocate_parallelepiped_domain

    subroutine destroy_and_deallocate_external_field(external_field)
        class(Abstract_External_Field), allocatable, intent(inout) :: external_field

        call external_field%destroy()
        if (allocated(external_field)) deallocate(external_field)
    end subroutine destroy_and_deallocate_external_field

    subroutine destroy_and_deallocate_reciprocal_lattice(reciprocal_lattice)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(inout) :: reciprocal_lattice

        call reciprocal_lattice%destroy()
        if (allocated(reciprocal_lattice)) deallocate(reciprocal_lattice)
    end subroutine destroy_and_deallocate_reciprocal_lattice

    subroutine deallocate_floor_penetration(floor_penetration)
        class(Abstract_Floor_Penetration), allocatable, intent(inout) :: floor_penetration

        if (allocated(floor_penetration)) deallocate(floor_penetration)
    end subroutine deallocate_floor_penetration

    subroutine destroy_and_deallocate_walls_potential(walls_potential)
        class(Abstract_Walls_Potential), allocatable, intent(inout) :: walls_potential

        call walls_potential%destroy()
        if (allocated(walls_potential)) deallocate(walls_potential)
    end subroutine destroy_and_deallocate_walls_potential

end module procedures_environment_factory
