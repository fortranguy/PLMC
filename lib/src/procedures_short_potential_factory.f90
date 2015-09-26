module procedures_short_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box, &
    XYZ_Periodic_Box, XY_Periodic_Box
use class_particles_diameter, only: Abstract_Particles_Diameter, &
    Null_Particles_Diameter
use class_particles_positions, only: Abstract_Particles_Positions
use types_particles, only: Particles_Wrapper
use procedures_types_selectors, only: particles_exist, particles_have_positions
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
use types_short_potential, only: Short_Potential_Wrapper

implicit none

private
public :: potential_factory_create, potential_factory_destroy, &
    allocate_and_set_expression, allocate_and_construct_pair, &
    allocate_and_construct_particles, allocate_list, allocate_and_construct_cells

contains

    subroutine potential_factory_create(short_potential, input_data, prefix, periodic_box, &
        particles)
        type(Short_Potential_Wrapper), intent(out) :: short_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Particles_Wrapper), intent(in) :: particles

        call allocate_and_set_expression(short_potential%expression, input_data, prefix, &
            particles%diameter)
        call allocate_and_construct_pair(short_potential%pair, input_data, prefix, &
            particles%diameter, short_potential%expression)
        call allocate_and_construct_particles(short_potential%particles, periodic_box, &
            particles%positions, short_potential%pair)
        call allocate_list(short_potential%list, input_data, prefix, particles%positions)
        call allocate_and_construct_cells(short_potential%cells, short_potential%list, &
            periodic_box, particles%positions, short_potential%pair)
    end subroutine potential_factory_create

    subroutine allocate_and_set_expression(potential_expression, input_data, prefix, &
        particles_diameter)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        call allocate_expression(potential_expression, input_data, prefix, particles_diameter)
        call set_expression(potential_expression, input_data, prefix)
    end subroutine allocate_and_set_expression

    subroutine allocate_expression(potential_expression, input_data, prefix, particles_diameter)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        if (particles_exist(particles_diameter)) then
            data_field = prefix//".Potential.name"
            call input_data%get(data_field, potential_name, data_found)
            call test_data_found(data_field, data_found)
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
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: LJ_epsilon, LJ_sigma

        select type(potential_expression)
            type is (Null_Potential_Expression)
                call potential_expression%set()
            type is (Lennard_Jones_Expression)
                data_field = prefix//".Potential.epsilon"
                call input_data%get(data_field, LJ_epsilon, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//".Potential.sigma"
                call input_data%get(data_field, LJ_sigma, data_found)
                call test_data_found(data_field, data_found)
                call potential_expression%set(LJ_epsilon, LJ_sigma)
            class default
                call error_exit("potential_expression type unknown.")
        end select
        if (allocated(data_field)) deallocate(data_field)
    end subroutine set_expression

    subroutine allocate_and_construct_pair(pair_potential, input_data, prefix, &
        particles_diameter, potential_expression)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression

        call allocate_pair(pair_potential, input_data, prefix, particles_diameter)
        call construct_pair(pair_potential, input_data, prefix, particles_diameter, &
            potential_expression)
    end subroutine allocate_and_construct_pair

    subroutine allocate_pair(pair_potential, input_data, prefix, particles_diameter)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (particles_exist(particles_diameter)) then
            data_field = prefix//".Potential.tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call test_data_found(data_field, data_found)
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

    subroutine construct_pair(pair_potential, input_data, prefix, particles_diameter, &
        potential_expression)
        class(Abstract_Pair_Potential), intent(inout) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain

        if (particles_exist(particles_diameter)) then
            domain%min = particles_diameter%get_min()
            data_field = prefix//".Potential.max distance"
            call input_data%get(data_field, domain%max, data_found)
            call test_data_found(data_field, data_found)
            select type (pair_potential)
                type is (Tabulated_Pair_Potential)
                    data_field = prefix//".Potential.delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call test_data_found(data_field, data_found)
            end select
        end if
        call pair_potential%construct(domain, potential_expression)
    end subroutine construct_pair

    subroutine allocate_and_construct_particles(particles_potential, periodic_box, &
        particles_positions, pair_potential)
        class(Abstract_Particles_Potential), allocatable, intent(out) :: particles_potential
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Positions), target, intent(in) :: particles_positions
        class(Abstract_Pair_Potential), target, intent(in) :: pair_potential

        if (particles_have_positions(particles_positions)) then
            allocate(Concrete_Particles_Potential :: particles_potential)
        else
            allocate(Null_Particles_Potential :: particles_potential)
        end if
        call particles_potential%construct(periodic_box, particles_positions, pair_potential)
    end subroutine allocate_and_construct_particles

    subroutine allocate_list(visitable_list, input_data, prefix, particles_positions)
        class(Abstract_Visitable_List), allocatable, intent(out) :: visitable_list
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Positions), intent(in) :: particles_positions

        character(len=:), allocatable :: data_field, cells_data_structure
        logical :: data_found

        if (particles_have_positions(particles_positions)) then
            data_field = prefix//".Potential.Cells.data structure"
            call input_data%get(data_field, cells_data_structure, data_found)
            call test_data_found(data_field, data_found)
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

    subroutine allocate_and_construct_cells(visitable_cells, visitable_list, periodic_box, &
        particles_positions, pair_potential)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: visitable_cells
        class(Abstract_Visitable_List), intent(in) :: visitable_list
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        if (particles_have_positions(particles_positions)) then
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

    subroutine potential_factory_destroy(short_potential)
        type(Short_Potential_Wrapper), intent(inout) :: short_potential

        call short_potential%cells%destroy()
        if (allocated(short_potential%cells)) deallocate(short_potential%cells)
        if (allocated(short_potential%list)) deallocate(short_potential%list)
        call short_potential%particles%destroy()
        if (allocated(short_potential%particles)) deallocate(short_potential%particles)
        call short_potential%pair%destroy()
        if (allocated(short_potential%pair)) deallocate(short_potential%pair)
        if (allocated(short_potential%expression)) deallocate(short_potential%expression)
    end subroutine potential_factory_destroy

end module procedures_short_potential_factory
