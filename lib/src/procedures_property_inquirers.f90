module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_walls_potential, only: Abstract_Walls_Potential, Null_Walls_Potential
use class_component_number, only: Abstract_Component_Number, Null_Component_Number
use class_component_diameter, only: Abstract_Component_Diameter, Null_Component_Diameter
use class_component_moment_norm, only: Abstract_Component_Moment_Norm, Null_Component_Moment_Norm
use class_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Null_Component_Dipolar_Moments
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Null_Component_Chemical_Potential
use class_moved_positions, only: Abstract_Moved_Positions, Null_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations, Null_Rotated_Orientations
use class_component_exchange, only: Abstract_Component_Exchange, Null_Component_Exchange
use class_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, Null_Ewald_Real_Pair

implicit none

private
public :: apply_external_field, use_reciprocal_lattice, use_walls, &
    component_exists, component_is_dipolar, component_has_positions, component_can_move, &
    component_has_orientations, component_can_exchange, component_can_rotate, component_interacts, &
    component_can_change

interface apply_external_field
    module procedure :: apply_external_field_from_json
end interface apply_external_field

interface use_reciprocal_lattice
    module procedure :: use_reciprocal_lattice_from_json
end interface use_reciprocal_lattice

interface use_walls
    module procedure :: use_walls_from_json
    module procedure :: use_walls_from_walls_potential
end interface use_walls

interface component_exists
    module procedure :: component_exists_from_json
    module procedure :: component_exists_from_number
end interface component_exists

interface component_is_dipolar
    module procedure :: component_is_dipolar_from_json
    module procedure :: component_is_dipolar_from_dipolar_moments
end interface component_is_dipolar

interface component_can_exchange
    module procedure :: component_can_exchange_from_json
    module procedure :: component_can_exchange_from_chemical_potential
    module procedure :: component_can_exchange_from_component_exchange
end interface component_can_exchange

interface component_interacts
    module procedure :: component_interacts_short
    module procedure :: component_interacts_long
end interface component_interacts

contains

    logical function apply_external_field_from_json(input_data, prefix) result(apply_external_field)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        apply_external_field = logical_from_json(input_data, prefix//"External Field.apply")
    end function apply_external_field_from_json

    logical function use_reciprocal_lattice_from_json(input_data, prefix) &
        result(use_reciprocal_lattice)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_reciprocal_lattice = logical_from_json(input_data, prefix//"Reciprocal Lattice.use")
    end function use_reciprocal_lattice_from_json

    logical function use_walls_from_json(input_data, prefix) result(use_walls)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_walls = logical_from_json(input_data, prefix//"Walls.use")
    end function use_walls_from_json

    logical function use_walls_from_walls_potential(walls_potential) result(use_walls)
        class(Abstract_Walls_Potential), intent(in) :: walls_potential

        select type (walls_potential)
            type is (Null_Walls_Potential)
                use_walls = .false.
            class default
                use_walls = .true.
        end select
    end function use_walls_from_walls_potential

    logical function component_exists_from_json(input_data, prefix) result(component_exists)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        component_exists = logical_from_json(input_data, prefix//"exists")
    end function component_exists_from_json

    pure logical function component_exists_from_number(component_number) result(component_exists)
        class(Abstract_Component_Number), intent(in) :: component_number

        select type (component_number)
            type is (Null_Component_Number)
                component_exists = .false.
            class default
                component_exists = .true.
        end select
    end function component_exists_from_number

    pure logical function component_interacts_short(pair_potential) result(component_interacts)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        select type (pair_potential)
            type is (Null_Pair_Potential)
                component_interacts = .false.
            class default
                component_interacts = .true.
        end select
    end function component_interacts_short

    pure logical function component_interacts_long(ewald_pair) result(component_interacts)
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_pair

        select type (ewald_pair)
            type is (Null_Ewald_Real_Pair)
                component_interacts = .false.
            class default
                component_interacts = .true.
        end select
    end function component_interacts_long

    pure logical function component_has_positions(partcles_positions)
        class(Abstract_Component_Coordinates), intent(in) :: partcles_positions

        select type (partcles_positions)
            type is (Concrete_Component_Positions)
                component_has_positions = .true.
            class default
                component_has_positions = .false.
        end select
    end function component_has_positions

    pure logical function component_can_move(moved_positions)
        class(Abstract_Moved_Positions), intent(in) :: moved_positions

        select type (moved_positions)
            type is (Null_Moved_Positions)
                component_can_move = .false.
            class default
                component_can_move = .true.
        end select
    end function component_can_move

    pure logical function component_can_rotate(rotated_orientations)
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations

        select type (rotated_orientations)
            type is (Null_Rotated_Orientations)
                component_can_rotate = .false.
            class default
                component_can_rotate = .true.
        end select
    end function component_can_rotate

    logical function component_is_dipolar_from_json(input_data, prefix) &
        result(component_is_dipolar)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (component_exists(input_data, prefix)) then
            component_is_dipolar = logical_from_json(input_data, prefix//"is dipolar")
        else
            component_is_dipolar = .false.
        end if
    end function component_is_dipolar_from_json

    pure logical function component_is_dipolar_from_moment_norm(component_moment_norm) &
        result(component_is_dipolar)
        class(Abstract_Component_Moment_Norm), intent(in) :: component_moment_norm

        select type (component_moment_norm)
            type is (Null_Component_Moment_Norm)
                component_is_dipolar = .false.
            class default
                component_is_dipolar = .true.
        end select
    end function component_is_dipolar_from_moment_norm

    pure logical function component_has_orientations(component_orientations)
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations

        select type (component_orientations)
            type is (Concrete_Component_Orientations)
                component_has_orientations = .true.
            class default
                component_has_orientations = .false.
        end select
    end function component_has_orientations

    pure logical function component_is_dipolar_from_dipolar_moments(component_dipolar_moments) &
        result(component_is_dipolar)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments

        select type (component_dipolar_moments)
            type is (Null_Component_Dipolar_Moments)
                component_is_dipolar = .false.
            class default
                component_is_dipolar = .true.
        end select
    end function component_is_dipolar_from_dipolar_moments

    logical function component_can_exchange_from_json(input_data, prefix) &
        result(component_can_exchange)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (component_exists(input_data, prefix)) then
            component_can_exchange = logical_from_json(input_data, prefix//"can exchange")
        else
            component_can_exchange = .false.
        end if
    end function component_can_exchange_from_json

    pure logical function component_can_exchange_from_chemical_potential(&
        component_chemical_potential) result(component_can_exchange)
        class(Abstract_Component_Chemical_Potential), intent(in) :: component_chemical_potential

        select type (component_chemical_potential)
            type is (Null_Component_Chemical_Potential)
                component_can_exchange = .false.
            class default
                component_can_exchange = .true.
        end select
    end function component_can_exchange_from_chemical_potential

    pure logical function component_can_exchange_from_component_exchange(component_exchange) &
        result(component_can_exchange)
        class(Abstract_Component_Exchange), intent(in) :: component_exchange

        select type (component_exchange)
            type is (Null_Component_Exchange)
                component_can_exchange = .false.
            class default
                component_can_exchange = .true.
        end select
    end function component_can_exchange_from_component_exchange

    pure logical function component_can_change(moved_positions, rotated_orientations, &
        component_exchange)
        class(Abstract_Moved_Positions), intent(in) :: moved_positions
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations
        class(Abstract_Component_Exchange), intent(in) :: component_exchange

        component_can_change = component_can_move(moved_positions) .or. &
            component_can_rotate(rotated_orientations) .or. &
            component_can_exchange(component_exchange)
    end function component_can_change

    logical function logical_from_json(input_data, switch_name)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: switch_name

        logical :: data_found

        call input_data%get(switch_name, logical_from_json, data_found)
        call check_data_found(switch_name, data_found)
    end function logical_from_json

end module procedures_property_inquirers
