module procedures_short_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_wrappers_prefix, only: environment_prefix
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box, &
    XYZ_Periodic_Box, XY_Periodic_Box
use class_floor_penetration, only: Abstract_Floor_Penetration
use types_environment_wrapper, only: Environment_Wrapper
use class_particles_diameter, only: Abstract_Particles_Diameter, &
    Null_Particles_Diameter
use class_component_positions, only: Abstract_Component_Positions
use types_particles_wrapper, only: Particles_Wrapper
use class_potential_expression, only: Abstract_Potential_Expression, &
    Null_Potential_Expression, Lennard_Jones_Expression
use types_potential_domain, only: Concrete_Potential_Domain
use class_pair_potential, only: Abstract_Pair_Potential, &
    Null_Pair_Potential, Tabulated_Pair_Potential, Raw_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential, &
    Concrete_Particles_Potential, Null_Particles_Potential
use class_visitable_list, only: Abstract_Visitable_List, &
    Null_Visitable_List, Concrete_Visitable_List, Concrete_Visitable_Array
use class_visitable_cells, only: Abstract_Visitable_Cells, &
    Null_Visitable_Cells, XYZ_PBC_Visitable_Cells, XY_PBC_Visitable_Cells
use types_short_potential_wrapper, only: Short_Potential_Wrapper, Short_Potential_Macro_Wrapper
use procedures_property_inquirers, only: use_walls, particles_exist, particles_have_positions, &
    particles_interact

implicit none

private
public :: short_potential_factory_create, short_potential_factory_destroy

interface short_potential_factory_create
    module procedure :: short_potential_factory_create_all
    module procedure :: short_potential_factory_create_macro
    module procedure :: short_potential_factory_create_inter_pair
    module procedure :: allocate_and_set_expression
    module procedure :: allocate_and_construct_pair
    module procedure :: allocate_and_construct_particles
    module procedure :: allocate_list
    module procedure :: allocate_and_construct_cells
end interface short_potential_factory_create

interface short_potential_factory_destroy
    module procedure :: destroy_and_deallocate_cells
    module procedure :: deallocate_list
    module procedure :: destroy_and_deallocate_particles
    module procedure :: destroy_and_deallocate_pair
    module procedure :: deallocate_expression
    module procedure :: short_potential_factory_destroy_macro
    module procedure :: short_potential_factory_destroy_all
end interface short_potential_factory_destroy

contains

    subroutine short_potential_factory_create_all(short_potential, environment, particles, &
        input_data, prefix)
        type(Short_Potential_Wrapper), intent(out) :: short_potential
        type(Environment_Wrapper), intent(in) :: environment
        type(Particles_Wrapper), intent(in) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Potential_Expression), allocatable :: expression, wall_expression
        class(Abstract_Visitable_List), allocatable :: list
        logical :: exist, exist_and_walls_used, interact

        exist = particles_exist(particles%number)
        call short_potential_factory_create(expression, exist, input_data, prefix)
        call short_potential_factory_create(short_potential%pair, exist, particles%diameter, &
            expression, input_data, prefix)
        call short_potential_factory_destroy(expression)
        exist_and_walls_used = exist .and. use_walls(input_data, environment_prefix)
        call short_potential_factory_create(wall_expression, exist_and_walls_used, &
            input_data, prefix//"With Walls.")
        call short_potential_factory_create(short_potential%wall_pair, exist_and_walls_used, &
            particles%wall_diameter, wall_expression, input_data, prefix//"With Walls.")
        call short_potential_factory_destroy(wall_expression)
        call short_potential_factory_create(short_potential%particles, environment%periodic_box, &
            particles%positions)
        interact = particles_interact(short_potential%pair)
        call short_potential_factory_create(list, interact, input_data, prefix)
        call short_potential_factory_create(short_potential%cells, interact, list, &
            environment%periodic_box, particles%positions, short_potential%pair)
        call short_potential_factory_destroy(list)
    end subroutine short_potential_factory_create_all

    subroutine short_potential_factory_create_macro(short_potential_macro, inter_pair, &
        periodic_box, particles_positions, input_data, prefix)
        type(Short_Potential_Macro_Wrapper), intent(out) :: short_potential_macro
        class(Abstract_Pair_Potential), intent(in) :: inter_pair
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Positions), intent(in) :: particles_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Visitable_List), allocatable :: list
        logical :: interact

        interact = particles_interact(inter_pair)
        call short_potential_factory_create(list, interact, input_data, prefix)
        call short_potential_factory_create(short_potential_macro%cells, interact, list, &
            periodic_box, particles_positions, inter_pair)
        call short_potential_factory_destroy(list)
    end subroutine short_potential_factory_create_macro

    subroutine short_potential_factory_create_inter_pair(inter_pair, exist, &
        particles_diameter, input_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: inter_pair
        logical, intent(in) :: exist
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Potential_Expression), allocatable :: expression

        call short_potential_factory_create(expression, exist, input_data, prefix)
        call short_potential_factory_create(inter_pair, exist, particles_diameter, &
            expression, input_data, prefix)
        call short_potential_factory_destroy(expression)
    end subroutine short_potential_factory_create_inter_pair

    subroutine allocate_and_set_expression(potential_expression, exist, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_expression(potential_expression, exist, input_data, prefix)
        call set_expression(potential_expression, input_data, prefix)
    end subroutine allocate_and_set_expression

    subroutine allocate_expression(potential_expression, exist, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        if (exist) then
            data_field = prefix//"name"
            call input_data%get(data_field, potential_name, data_found)
            call check_data_found(data_field, data_found)
            select case(potential_name)
                case ("null")
                    allocate(Null_Potential_Expression :: potential_expression)
                case ("LJ")
                    allocate(Lennard_Jones_Expression :: potential_expression)
                case default
                    call error_exit(potential_name//" unknown potential_name."//&
                        "Choose between: 'null' and LJ.")
            end select
            deallocate(potential_name)
            deallocate(data_field)
        else
            allocate(Null_Potential_Expression :: potential_expression)
        end if
    end subroutine allocate_expression

    subroutine set_expression(potential_expression, input_data, prefix)
        class(Abstract_Potential_Expression), intent(inout) :: potential_expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: LJ_epsilon, LJ_sigma

        select type(potential_expression)
            type is (Null_Potential_Expression)
                call potential_expression%set()
            type is (Lennard_Jones_Expression)
                data_field = prefix//"epsilon"
                call input_data%get(data_field, LJ_epsilon, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"sigma"
                call input_data%get(data_field, LJ_sigma, data_found)
                call check_data_found(data_field, data_found)
                call potential_expression%set(LJ_epsilon, LJ_sigma)
            class default
                call error_exit("potential_expression type unknown.")
        end select
        if (allocated(data_field)) deallocate(data_field)
    end subroutine set_expression

    subroutine allocate_and_construct_pair(pair_potential, exist, particles_diameter, &
        potential_expression, input_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        logical, intent(in) :: exist
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_pair(pair_potential, exist, input_data, prefix)
        call construct_pair(pair_potential, exist, particles_diameter, potential_expression, &
            input_data, prefix)
    end subroutine allocate_and_construct_pair

    subroutine allocate_pair(pair_potential, exist, input_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (exist) then
            data_field = prefix//"tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_Pair_Potential :: pair_potential)
            else
                allocate(Raw_Pair_Potential :: pair_potential)
            end if
            deallocate(data_field)
        else
            allocate(Null_Pair_Potential :: pair_potential)
        end if
    end subroutine allocate_pair

    subroutine construct_pair(pair_potential, exist, particles_diameter, potential_expression, &
        input_data, prefix)
        class(Abstract_Pair_Potential), intent(inout) :: pair_potential
        logical, intent(in) :: exist
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain

        if (exist) then
            domain%min = particles_diameter%get_min()
            select type (potential_expression)
                type is (Null_Potential_Expression)
                    domain%max = domain%min
                class default
                    data_field = prefix//"maximum distance"
                    call input_data%get(data_field, domain%max, data_found)
                    call check_data_found(data_field, data_found)
            end select
            select type (pair_potential)
                type is (Tabulated_Pair_Potential)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair_potential%construct(domain, potential_expression)
    end subroutine construct_pair

    subroutine allocate_and_construct_particles(particles_potential, periodic_box, &
        particles_positions)
        class(Abstract_Particles_Potential), allocatable, intent(out) :: particles_potential
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Positions), intent(in) :: particles_positions

        if (particles_have_positions(particles_positions)) then
            allocate(Concrete_Particles_Potential :: particles_potential)
        else
            allocate(Null_Particles_Potential :: particles_potential)
        end if
        call particles_potential%construct(periodic_box, particles_positions)
    end subroutine allocate_and_construct_particles

    subroutine allocate_list(visitable_list, interact, input_data, prefix)
        class(Abstract_Visitable_List), allocatable, intent(out) :: visitable_list
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, cells_data_structure
        logical :: data_found

        if (interact) then
            data_field = prefix//"Cells.data structure"
            call input_data%get(data_field, cells_data_structure, data_found)
            call check_data_found(data_field, data_found)
            select case(cells_data_structure)
                case ("list")
                    allocate(Concrete_Visitable_List :: visitable_list)
                case ("array")
                    allocate(Concrete_Visitable_Array :: visitable_list)
                case default
                    call error_exit(cells_data_structure//" unknown."&
                        //"Choose between 'list' and 'array'.")
            end select
            deallocate(cells_data_structure)
        else
            allocate(Null_Visitable_List :: visitable_list)
        end if
    end subroutine allocate_list

    subroutine allocate_and_construct_cells(visitable_cells, interact, visitable_list, &
        periodic_box, particles_positions, pair_potential)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells
        logical, intent(in) :: interact
        class(Abstract_Visitable_List), intent(in) :: visitable_list
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Positions), intent(in) :: particles_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        if (interact) then
            select type(periodic_box)
                type is (XYZ_Periodic_Box)
                    allocate(XYZ_PBC_Visitable_Cells :: visitable_cells)
                type is (XY_Periodic_Box)
                    allocate(XY_PBC_Visitable_Cells :: visitable_cells)
                class default
                    call error_exit("periodic_box type unknown.")
            end select
        else
            allocate(Null_Visitable_Cells :: visitable_cells)
        end if
        call visitable_cells%construct(visitable_list, periodic_box, particles_positions, &
            pair_potential)
    end subroutine allocate_and_construct_cells

    subroutine short_potential_factory_destroy_macro(short_potential_macro)
        type(Short_Potential_Macro_Wrapper), intent(inout) :: short_potential_macro

        call short_potential_factory_destroy(short_potential_macro%cells)
    end subroutine short_potential_factory_destroy_macro

    subroutine short_potential_factory_destroy_all(short_potential)
        type(Short_Potential_Wrapper), intent(inout) :: short_potential

        call short_potential_factory_destroy(short_potential%cells)
        call short_potential_factory_destroy(short_potential%particles)
        call short_potential_factory_destroy(short_potential%wall_pair)
        call short_potential_factory_destroy(short_potential%pair)
    end subroutine short_potential_factory_destroy_all

    subroutine deallocate_expression(potential_expression)
        class(Abstract_Potential_Expression), allocatable, intent(inout) :: potential_expression

        if (allocated(potential_expression)) deallocate(potential_expression)
    end subroutine deallocate_expression

    subroutine destroy_and_deallocate_pair(pair_potential)
        class(Abstract_Pair_Potential), allocatable, intent(inout) :: pair_potential

        call pair_potential%destroy()
        if (allocated(pair_potential)) deallocate(pair_potential)
    end subroutine destroy_and_deallocate_pair

    subroutine destroy_and_deallocate_particles(particles_potential)
        class(Abstract_Particles_Potential), allocatable, intent(inout) :: particles_potential

        call particles_potential%destroy()
        if (allocated(particles_potential)) deallocate(particles_potential)
    end subroutine destroy_and_deallocate_particles

    subroutine deallocate_list(visitable_list)
        class(Abstract_Visitable_List), allocatable, intent(inout) :: visitable_list

        if (allocated(visitable_list)) deallocate(visitable_list)
    end subroutine deallocate_list

    subroutine destroy_and_deallocate_cells(visitable_cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: visitable_cells

        call visitable_cells%destroy()
        if (allocated(visitable_cells)) deallocate(visitable_cells)
    end subroutine destroy_and_deallocate_cells

end module procedures_short_potential_factory
