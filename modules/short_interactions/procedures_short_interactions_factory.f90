module procedures_short_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_number_to_string, only: Concrete_Number_to_String
use class_periodic_box, only: Abstract_Periodic_Box, &
    XYZ_Periodic_Box, XY_Periodic_Box
use class_walls_potential, only: Abstract_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use class_minimum_distance, only: Abstract_Minimum_Distance
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Minimum_Distance_Wrapper, Minimum_Distances_Wrapper, &
    Mixture_Wrapper
use class_potential_expression, only: Abstract_Potential_Expression, &
    Null_Potential_Expression, Lennard_Jones_Expression
use types_potential_domain, only: Short_Potential_Domain
use class_pair_potential, only: Abstract_Pair_Potential, &
    Null_Pair_Potential, Tabulated_Pair_Potential, Raw_Pair_Potential
use class_visitable_list, only: Abstract_Visitable_List, &
    Null_Visitable_List, Concrete_Visitable_List, Concrete_Visitable_Array
use class_visitable_cells, only: Abstract_Visitable_Cells, &
    Null_Visitable_Cells, XYZ_PBC_Visitable_Cells, XY_PBC_Visitable_Cells
use class_walls_potential_visitor, only: Abstract_Walls_Potential_Visitor, &
    Concrete_Walls_Potential_Visitor, Null_Walls_Potential_Visitor
use class_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor, &
    Concrete_Short_Pairs_Visitor, Null_Short_Pairs_Visitor
use types_short_interactions_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Wrapper, &
    Short_Interactions_Wrapper
use procedures_property_inquirers, only: use_walls, components_interact

implicit none

private
public :: short_interactions_create, short_interactions_destroy

interface short_interactions_create
    module procedure :: create_all
    module procedure :: create_walls_visitor
    module procedure :: create_components_visitor
    module procedure :: create_components_cells
    module procedure :: create_components_pairs
    module procedure :: create_wall_pairs
    module procedure :: create_pair
    module procedure :: create_expression
end interface short_interactions_create

interface short_interactions_destroy
    module procedure :: destroy_expression
    module procedure :: destroy_pair
    module procedure :: destroy_pairs
    module procedure :: destroy_components_pairs
    module procedure :: destroy_components_cells
    module procedure :: destroy_components_visitor
    module procedure :: destroy_walls_visitor
    module procedure :: destroy_all
end interface short_interactions_destroy

contains

    subroutine create_all(short_interactions, environment, mixture, input_data, prefix)
        type(Short_Interactions_Wrapper), intent(out) :: short_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: interact, interact_with_walls
        class(Abstract_Visitable_List), allocatable :: list

        call short_interactions_create(short_interactions%wall_pairs, interact_with_walls, &
            mixture%wall_min_distances, input_data, prefix)
        call short_interactions_create(short_interactions%walls_visitor, environment%&
            walls_potential, interact_with_walls)
        call short_interactions_create(short_interactions%components_pairs, interact, mixture%&
            components_min_distances, input_data, prefix)
        call short_interactions_create(short_interactions%components_visitor, environment%&
            periodic_box, interact)
        call allocate_list(list, interact, input_data, prefix)
        call short_interactions_create(short_interactions%components_cells, environment%&
            periodic_box, mixture%components, short_interactions%components_pairs, interact, list)
        call deallocate_list(list)
    end subroutine create_all

    subroutine destroy_all(short_interactions)
        type(Short_Interactions_Wrapper), intent(inout) :: short_interactions

        call short_interactions_destroy(short_interactions%walls_visitor)
        call short_interactions_destroy(short_interactions%components_visitor)
        call short_interactions_destroy(short_interactions%components_cells)
        call short_interactions_destroy(short_interactions%wall_pairs)
        call short_interactions_destroy(short_interactions%components_pairs)
    end subroutine destroy_all

    subroutine create_walls_visitor(visitor, walls_potential, interact)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Walls_Potential), intent(in) :: walls_potential
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Walls_Potential_Visitor :: visitor)
        else
            allocate(Null_Walls_Potential_Visitor :: visitor)
        end if
        call visitor%construct(walls_potential)
    end subroutine create_walls_visitor

    subroutine destroy_walls_visitor(visitor)
        class(Abstract_Walls_Potential_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_walls_visitor

    subroutine create_wall_pairs(pairs, interact, min_distances, input_data, prefix)
        type(Pair_Potential_Wrapper), allocatable, intent(out) :: pairs(:)
        logical, intent(out) :: interact
        type(Minimum_Distance_Wrapper), intent(in) :: min_distances(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        class(Abstract_Potential_Expression), allocatable :: expression
        character(len=:), allocatable :: pair_prefix
        type(Concrete_Number_to_String) :: string

        interact = .true.
        allocate(pairs(size(min_distances)))
        do i_component = 1, size(pairs)
            pair_prefix = prefix//"Component "//string%get(i_component)//".With Walls."
            associate (min_distance => min_distances(i_component)%min_distance, &
                interact_i => components_interact(min_distances(i_component)%min_distance))
                interact = interact .and. interact_i
                call short_interactions_create(expression, interact_i, input_data, pair_prefix)
                call short_interactions_create(pairs(i_component)%pair_potential, min_distance, &
                    interact_i, expression, input_data, pair_prefix)
            end associate
            deallocate(pair_prefix)
            call short_interactions_destroy(expression)
        end do
    end subroutine create_wall_pairs

    subroutine create_components_visitor(visitor, periodic_box, interact)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: interact

        if (interact) then
            allocate(Concrete_Short_Pairs_Visitor :: visitor)
        else
            allocate(Null_Short_Pairs_Visitor :: visitor)
        end if
        call visitor%construct(periodic_box)
    end subroutine create_components_visitor

    subroutine destroy_components_visitor(visitor)
        class(Abstract_Short_Pairs_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_components_visitor

    subroutine create_components_cells(cells, periodic_box, components, pairs, interact, list)
        class(Abstract_Visitable_Cells), allocatable, intent(out) :: cells(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(Pair_Potentials_Wrapper), intent(in) :: pairs(:)
        logical, intent(in) :: interact
        class(Abstract_Visitable_List), intent(in) :: list

        integer :: j_component, i_component
        integer :: j_pair, i_pair

        if (interact) then
            select type(periodic_box)
                type is (XYZ_Periodic_Box)
                    allocate(XYZ_PBC_Visitable_Cells :: cells(size(pairs), size(pairs)))
                type is (XY_Periodic_Box)
                    allocate(XY_PBC_Visitable_Cells :: cells(size(pairs), size(pairs)))
                class default
                    call error_exit("periodic_box type unknown.")
            end select
        else
            allocate(Null_Visitable_Cells :: cells(size(pairs), size(pairs)))
        end if

        do j_component = 1, size(cells, 2)
            do i_component = 1, size(cells, 1)
                j_pair = maxval([j_component, i_component])
                i_pair = minval([j_component, i_component])
                associate (pair_ij => pairs(j_pair)%with_components(i_pair)%pair_potential)
                    call cells(i_component, j_component)%construct(list, periodic_box, &
                        components(i_component)%positions, pair_ij)
                end associate
            end do
        end do
    end subroutine create_components_cells

    subroutine destroy_components_cells(cells)
        class(Abstract_Visitable_Cells), allocatable, intent(inout) :: cells(:, :)

        integer :: j_component, i_component

        if (allocated(cells)) then
            do j_component = size(cells, 2), 1, -1
                do i_component = size(cells, 1), 1, -1
                    call cells(i_component, j_component)%destroy()
                end do
            end do
            deallocate(cells)
        end if
    end subroutine destroy_components_cells

    subroutine allocate_list(list, interact, input_data, prefix)
        class(Abstract_Visitable_List), allocatable, intent(out) :: list
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
                    allocate(Concrete_Visitable_List :: list)
                case ("array")
                    allocate(Concrete_Visitable_Array :: list)
                case default
                    call error_exit(cells_data_structure//" unknown."&
                        //"Choose between 'list' and 'array'.")
            end select
            deallocate(cells_data_structure)
        else
            allocate(Null_Visitable_List :: list)
        end if
    end subroutine allocate_list

    subroutine deallocate_list(list)
        class(Abstract_Visitable_List), allocatable, intent(inout) :: list

        if (allocated(list)) deallocate(list)
    end subroutine deallocate_list

    subroutine create_components_pairs(pairs, interact, min_distances, input_data, prefix)
        type(Pair_Potentials_Wrapper), allocatable, intent(out) :: pairs(:)
        logical, intent(out) :: interact
        type(Minimum_Distances_Wrapper), intent(in) :: min_distances(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: j_component, i_component
        logical :: interact_ij
        class(Abstract_Potential_Expression), allocatable :: expression
        character(len=:), allocatable :: pair_prefix
        type(Concrete_Number_to_String) :: string

        interact = .true.
        allocate(pairs(size(min_distances)))
        do j_component = 1, size(pairs)
            allocate(pairs(j_component)%with_components(j_component))
            do i_component = 1, size(pairs(j_component)%with_components)
                if (i_component == j_component) then
                    pair_prefix = prefix//"Component "//string%get(i_component)//"."
                else
                    pair_prefix = prefix//"Inter "//string%get(i_component)//&
                        string%get(j_component)//"."
                end if
                associate (min_distance => min_distances(j_component)%with_components(i_component)%&
                        min_distance)
                    interact_ij = components_interact(min_distance)
                    interact = interact .and. interact_ij
                    call short_interactions_create(expression, interact_ij, input_data, &
                        pair_prefix)
                    call short_interactions_create(pairs(j_component)%&
                        with_components(i_component)%pair_potential, min_distance, interact_ij, &
                        expression, input_data, pair_prefix)
                end associate
                deallocate(pair_prefix)
                call short_interactions_destroy(expression)
            end do
        end do
    end subroutine create_components_pairs

    subroutine destroy_components_pairs(pairs)
        type(Pair_Potentials_Wrapper), allocatable, intent(inout) :: pairs(:)

        integer :: i_component

        if (allocated(pairs)) then
            do i_component = size(pairs), 1, -1
                call short_interactions_destroy(pairs(i_component)%with_components)
            end do
            deallocate(pairs)
        end if
    end subroutine destroy_components_pairs

    subroutine destroy_pairs(pairs)
        type(Pair_Potential_Wrapper), allocatable, intent(inout) :: pairs(:)

        integer :: i_component

        if (allocated(pairs)) then
            do i_component = size(pairs), 1, -1
                call short_interactions_destroy(pairs(i_component)%pair_potential)
            end do
            deallocate(pairs)
        end if
    end subroutine destroy_pairs

    subroutine create_pair(pair, min_distance, interact, expression, input_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_Potential_Expression), intent(in) :: expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_pair(pair, interact, input_data, prefix)
        call construct_pair(pair, min_distance, interact, expression, input_data, prefix)
    end subroutine create_pair

    subroutine allocate_pair(pair, interact, input_data, prefix)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (interact) then
            data_field = prefix//"tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_Pair_Potential :: pair)
            else
                allocate(Raw_Pair_Potential :: pair)
            end if
        else
            allocate(Null_Pair_Potential :: pair)
        end if
    end subroutine allocate_pair

    subroutine construct_pair(pair, min_distance, interact, expression, input_data, prefix)
        class(Abstract_Pair_Potential), intent(inout) :: pair
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_Potential_Expression), intent(in) :: expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Short_Potential_Domain) :: domain

        if (interact) then
            domain%min = min_distance%get()
            select type (expression)
                type is (Null_Potential_Expression)
                    domain%max = domain%min
                class default
                    data_field = prefix//"maximum distance"
                    call input_data%get(data_field, domain%max, data_found)
                    call check_data_found(data_field, data_found)
            end select
            select type (pair)
                type is (Tabulated_Pair_Potential)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair%construct(domain, expression)
    end subroutine construct_pair

    subroutine destroy_pair(pair)
        class(Abstract_Pair_Potential), allocatable, intent(inout) :: pair

        call pair%destroy()
        if (allocated(pair)) deallocate(pair)
    end subroutine destroy_pair

    subroutine create_expression(expression, interact, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_expression(expression, interact, input_data, prefix)
        call set_expression(expression, input_data, prefix)
    end subroutine create_expression

    subroutine allocate_expression(expression, interact, input_data, prefix)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: expression
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        if (interact) then
            data_field = prefix//"name"
            call input_data%get(data_field, potential_name, data_found)
            call check_data_found(data_field, data_found)
            select case(potential_name)
                case ("null")
                    allocate(Null_Potential_Expression :: expression)
                case ("LJ")
                    allocate(Lennard_Jones_Expression :: expression)
                case default
                    call error_exit(potential_name//" unknown potential_name."//&
                        "Choose between: 'null' and LJ.")
            end select
        else
            allocate(Null_Potential_Expression :: expression)
        end if
    end subroutine allocate_expression

    subroutine set_expression(expression, input_data, prefix)
        class(Abstract_Potential_Expression), intent(inout) :: expression
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: LJ_epsilon, LJ_sigma

        select type(expression)
            type is (Null_Potential_Expression)
                call expression%set()
            type is (Lennard_Jones_Expression)
                data_field = prefix//"epsilon"
                call input_data%get(data_field, LJ_epsilon, data_found)
                call check_data_found(data_field, data_found)
                data_field = prefix//"sigma"
                call input_data%get(data_field, LJ_sigma, data_found)
                call check_data_found(data_field, data_found)
                call expression%set(LJ_epsilon, LJ_sigma)
            class default
                call error_exit("expression type unknown.")
        end select
    end subroutine set_expression

    subroutine destroy_expression(expression)
        class(Abstract_Potential_Expression), allocatable, intent(inout) :: expression

        if (allocated(expression)) deallocate(expression)
    end subroutine destroy_expression

end module procedures_short_interactions_factory
