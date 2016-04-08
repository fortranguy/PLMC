module procedures_environment_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use class_temperature, only: Abstract_Temperature, Concrete_Temperature
use class_field_expression, only: Abstract_Field_Expression, Constant_Field_Expression, &
    Centered_Plates_Expression, Null_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain, &
    Concrete_Parallelepiped_Domain, Concrete_Box_Domain, Null_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field, Concrete_External_Field, &
    Null_External_Field
use class_permittivity, only: Abstract_Permittivity, Concrete_Permittivity, Null_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice, &
    Null_Reciprocal_Lattice
use class_floor_penetration, only: Abstract_Floor_Penetration, Flat_Floor_Penetration, &
    Centered_Block_Penetration, Null_Floor_Penetration
use procedures_component_factory, only: component_destroy
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use procedures_short_interactions_factory, only: short_interactions_create, &
    short_interactions_destroy
use class_walls_potential, only: Abstract_Walls_Potential, Concrete_Walls_Potential, &
    Null_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use procedures_property_inquirers, only: periodicity_is_xyz, periodicity_is_xy, &
    apply_external_field, use_permittivity, use_reciprocal_lattice, use_walls

implicit none

private
public :: environment_create, environment_destroy

interface environment_create
    module procedure :: create_all
    module procedure :: create_periodic_box
    module procedure :: create_temperature
    module procedure :: create_field_expression
    module procedure :: create_parallelepiped_domain
    module procedure :: create_external_field
    module procedure :: create_permittivity
    module procedure :: create_reciprocal_lattice
    module procedure :: create_floor_penetration
    module procedure :: create_walls_potential
end interface environment_create

interface environment_destroy
    module procedure :: destroy_walls_potential
    module procedure :: destroy_floor_penetration
    module procedure :: destroy_reciprocal_lattice
    module procedure :: destroy_permittivity
    module procedure :: destroy_external_field
    module procedure :: destroy_parallelepiped_domain
    module procedure :: destroy_field_expression
    module procedure :: destroy_temperature
    module procedure :: destroy_periodic_box
    module procedure :: destroy_all
end interface environment_destroy

interface environment_check
    module procedure :: check_consistency
end interface environment_check

contains

    subroutine create_all(environment, input_data, prefix)
        type(Environment_Wrapper), intent(out) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        logical :: field_applied, walls_used

        call environment_create(environment%periodic_box, input_data, prefix)
        call environment_create(environment%temperature, input_data, prefix)
        call environment_create(environment%permittivity, input_data, prefix)
        field_applied = apply_external_field(input_data, prefix)
        call environment_create(field_expression, environment%permittivity, field_applied, &
            input_data, prefix)
        call environment_create(parallelepiped_domain, field_applied, &
            environment%periodic_box, input_data, prefix//"External Field.")
        call environment_create(environment%external_field, field_applied, &
            parallelepiped_domain, field_expression)
        call environment_destroy(field_expression)
        call environment_destroy(parallelepiped_domain)
        call environment_create(environment%reciprocal_lattice, environment%periodic_box, &
            input_data, prefix)
        walls_used = use_walls(input_data, prefix)
        call environment_create(floor_penetration, walls_used, input_data, prefix)
        call environment_create(environment%walls_potential, walls_used, &
            environment%periodic_box, floor_penetration, input_data, prefix)
        call environment_destroy(floor_penetration)

        call environment_check(environment%periodic_box, environment%walls_potential)
    end subroutine create_all

    subroutine destroy_all(environment)
        type(Environment_Wrapper), intent(inout) :: environment

        call environment_destroy(environment%walls_potential)
        call environment_destroy(environment%reciprocal_lattice)
        call environment_destroy(environment%external_field)
        call environment_destroy(environment%permittivity)
        call environment_destroy(environment%temperature)
        call environment_destroy(environment%periodic_box)
    end subroutine destroy_all

    subroutine create_periodic_box(periodic_box, input_data, prefix)
        class(Abstract_Periodic_Box), allocatable, intent(out) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: box_periodicity
        real(DP), allocatable :: box_size(:)
        character(len=:), allocatable :: data_field
        logical :: data_found

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
        call periodic_box%set(box_size)
    end subroutine create_periodic_box

    subroutine destroy_periodic_box(periodic_box)
        class(Abstract_Periodic_Box), allocatable, intent(inout) :: periodic_box

        if (allocated(periodic_box)) deallocate(periodic_box)
    end subroutine destroy_periodic_box

    subroutine create_temperature(temperature, input_data, prefix)
        class(Abstract_Temperature), allocatable, intent(out) :: temperature
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: temperature_value
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"Thermostat.temperature"
        call input_data%get(data_field, temperature_value, data_found)
        call check_data_found(data_field, data_found)
        allocate(Concrete_Temperature :: temperature)
        call temperature%set(temperature_value)
    end subroutine create_temperature

    subroutine destroy_temperature(temperature)
        class(Abstract_Temperature), allocatable, intent(inout) :: temperature

        if (allocated(temperature)) deallocate(temperature)
    end subroutine destroy_temperature

    subroutine create_field_expression(field_expression, permittivity, field_applied, input_data, &
        prefix)
        class(Abstract_Field_Expression), allocatable, intent(out) :: field_expression
        class(Abstract_Permittivity), intent(in) :: permittivity
        logical, intent(in) :: field_applied
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_field_expression(field_expression, field_applied, input_data, prefix)
        call set_field_expression(field_expression, permittivity, input_data, prefix)
    end subroutine create_field_expression

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

    subroutine destroy_field_expression(field_expression)
        class(Abstract_Field_Expression), allocatable, intent(inout) :: field_expression

        if (allocated(field_expression)) deallocate(field_expression)
    end subroutine destroy_field_expression

    subroutine create_parallelepiped_domain(parallelepiped_domain, needed,periodic_box, &
        input_data, prefix)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(out) :: parallelepiped_domain
        logical, intent(in) :: needed
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: domain_name, data_field
        real(DP), allocatable :: domain_origin(:), domain_size(:)
        logical :: data_found

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
                    allocate(Concrete_Box_Domain :: parallelepiped_domain)
                case default
                    call error_exit(domain_name//" domain_name unknown. "//&
                        "Choose between 'domain' and 'box'.")
            end select
        else
            allocate(Null_Parallelepiped_Domain :: parallelepiped_domain)
        end if
        call parallelepiped_domain%construct(periodic_box, domain_origin, domain_size)
        if (allocated(domain_size)) deallocate(domain_size)
        if (allocated(domain_origin)) deallocate(domain_origin)
    end subroutine create_parallelepiped_domain

    subroutine destroy_parallelepiped_domain(parallelepiped_domain)
        class(Abstract_Parallelepiped_Domain), allocatable, intent(inout) :: parallelepiped_domain

        call parallelepiped_domain%destroy()
        if (allocated(parallelepiped_domain)) deallocate(parallelepiped_domain)
    end subroutine destroy_parallelepiped_domain

    subroutine create_external_field(external_field, field_applied, parallelepiped_domain, &
        field_expression)
        class(Abstract_External_Field), allocatable, intent(out) :: external_field
        logical, intent(in) :: field_applied
        class(Abstract_Parallelepiped_Domain), intent(in) :: parallelepiped_domain
        class(Abstract_Field_Expression), intent(in) :: field_expression

        if (field_applied) then
            allocate(Concrete_External_Field :: external_field)
        else
            allocate(Null_External_Field :: external_field)
        end if
        call external_field%construct(parallelepiped_domain, field_expression)
    end subroutine create_external_field

    subroutine destroy_external_field(external_field)
        class(Abstract_External_Field), allocatable, intent(inout) :: external_field

        call external_field%destroy()
        if (allocated(external_field)) deallocate(external_field)
    end subroutine destroy_external_field

    subroutine create_permittivity(permittivity, input_data, prefix)
        class(Abstract_Permittivity), allocatable, intent(out) :: permittivity
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: permittivity_value

        if (use_permittivity(input_data, prefix)) then
            data_field = prefix//"Permittivity.value"
            call input_data%get(data_field, permittivity_value, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Permittivity :: permittivity)
        else
            allocate(Null_Permittivity :: permittivity)
        end if
        call permittivity%set(permittivity_value)
    end subroutine create_permittivity

    subroutine destroy_permittivity(permittivity)
        class(Abstract_Permittivity), allocatable, intent(inout) :: permittivity

        if (allocated(permittivity)) then
            deallocate(permittivity)
        end if
    end subroutine destroy_permittivity

    subroutine create_reciprocal_lattice(reciprocal_lattice, periodic_box, input_data, prefix)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(out) :: reciprocal_lattice
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer, allocatable :: numbers(:)

        if (use_reciprocal_lattice(input_data, prefix)) then
            data_field = prefix//"Reciprocal Lattice.numbers"
            call input_data%get(data_field, numbers, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Reciprocal_Lattice :: reciprocal_lattice)
        else
            allocate(Null_Reciprocal_Lattice :: reciprocal_lattice)
        end if
        call reciprocal_lattice%construct(periodic_box, numbers)
        if (allocated(numbers)) deallocate(numbers)
    end subroutine create_reciprocal_lattice

    subroutine destroy_reciprocal_lattice(reciprocal_lattice)
        class(Abstract_Reciprocal_Lattice), allocatable, intent(inout) :: reciprocal_lattice

        call reciprocal_lattice%destroy()
        if (allocated(reciprocal_lattice)) deallocate(reciprocal_lattice)
    end subroutine destroy_reciprocal_lattice

    subroutine create_floor_penetration(floor_penetration, walls_used, input_data, prefix)
        class(Abstract_Floor_Penetration), allocatable, intent(out) :: floor_penetration
        logical, intent(in) :: walls_used
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_floor_penetration(floor_penetration, walls_used, input_data, prefix)
        call set_floor_penetration(floor_penetration, input_data, prefix)
    end subroutine create_floor_penetration

    subroutine allocate_floor_penetration(floor_penetration, walls_used, input_data, prefix)
        class(Abstract_Floor_Penetration), allocatable, intent(out) :: floor_penetration
        logical, intent(in) :: walls_used
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: walls_name, data_field
        logical :: data_found

        if (walls_used) then
            data_field = prefix//"Walls.name"
            call input_data%get(data_field, walls_name, data_found)
            call check_data_found(data_field, data_found)
            select case(walls_name)
                case ("flat")
                    allocate(Flat_Floor_Penetration :: floor_penetration)
                case ("block")
                    allocate(Centered_Block_Penetration :: floor_penetration)
                case default
                    call error_exit(walls_name//" walls_name unknown. "//&
                        "Choose between 'flat' and 'block'.")
            end select
        else
            allocate(Null_Floor_Penetration :: floor_penetration)
        end if
    end subroutine allocate_floor_penetration

    subroutine set_floor_penetration(floor_penetration, input_data, prefix)
        class(Abstract_Floor_Penetration), intent(inout) :: floor_penetration
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP), allocatable :: block_size(:)
        real(DP) ::block_radius
        character(len=:), allocatable :: data_field
        logical :: data_found

        select type (floor_penetration)
            type is (Flat_Floor_Penetration)
            type is (Centered_Block_Penetration)
                data_field = prefix//"Walls.size"
                call input_data%get(data_field, block_size, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"Walls.radius"
                call input_data%get(data_field, block_radius, data_found)
                call check_data_found(data_field, data_found)
                call floor_penetration%set(block_size, block_radius)
            type is (Null_Floor_Penetration)
            class default
                call error_exit("floor_penetration type unknown.")
        end select
    end subroutine set_floor_penetration

    subroutine destroy_floor_penetration(floor_penetration)
        class(Abstract_Floor_Penetration), allocatable, intent(inout) :: floor_penetration

        if (allocated(floor_penetration)) deallocate(floor_penetration)
    end subroutine destroy_floor_penetration

    subroutine create_walls_potential(walls_potential, walls_used, periodic_box, &
        floor_penetration, input_data, prefix)
        class(Abstract_Walls_Potential), allocatable, intent(out) :: walls_potential
        logical, intent(in) :: walls_used
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: gap

        if (walls_used) then
            allocate(Concrete_Walls_Potential :: walls_potential)
            data_field = prefix//"Walls.gap"
            call input_data%get(data_field, gap, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Walls_Potential :: walls_potential)
        end if
        call walls_potential%construct(periodic_box, gap, floor_penetration)
    end subroutine create_walls_potential

    subroutine destroy_walls_potential(walls_potential)
        class(Abstract_Walls_Potential), allocatable, intent(inout) :: walls_potential

        call walls_potential%destroy()
        if (allocated(walls_potential)) deallocate(walls_potential)
    end subroutine destroy_walls_potential

    subroutine check_consistency(periodic_box, walls_potential)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Walls_Potential), intent(in) :: walls_potential

        if (periodicity_is_xyz(periodic_box) .and. use_walls(walls_potential)) then
            call warning_continue("environment_check: periodicity is XYZ but walls are used.")
        end if
        if (periodicity_is_xy(periodic_box) .and. .not.use_walls(walls_potential)) then
            call warning_continue("environment_check: periodicity is XY but walls are not used.")
        end if
    end subroutine check_consistency

end module procedures_environment_factory
