module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_permittivity, only: Abstract_Permittivity, Concrete_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice
use class_walls_potential, only: Abstract_Walls_Potential, Concrete_Walls_Potential
use class_component_number, only: Abstract_Component_Number, Concrete_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations
use class_minimum_distance, only: Abstract_Minimum_Distance, Concrete_Minimum_Distance
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Concrete_Component_Dipolar_Moments
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential
use class_changed_coordinates, only: Abstract_Changed_Coordinates
use class_moved_positions, only: Concrete_Moved_Positions
use class_rotated_orientations, only: Concrete_Rotated_Orientations
use class_component_exchange, only: Abstract_Component_Exchange, Concrete_Component_Exchange
use class_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, Null_Ewald_Real_Pair

implicit none

private
public :: apply_external_field, use_permittivity, use_reciprocal_lattice, use_walls, &
    component_exists, component_is_dipolar, component_has_positions, component_can_move, &
    component_has_orientations, component_can_exchange, component_can_rotate, components_interact, &
    component_can_change

interface apply_external_field
    module procedure :: apply_external_field_from_json
end interface apply_external_field

interface use_permittivity
    module procedure :: use_permittivity_from_json
    module procedure :: use_permittivity_from
end interface

interface use_reciprocal_lattice
    module procedure :: use_reciprocal_lattice_from_json
    module procedure :: use_reciprocal_lattice_from
end interface use_reciprocal_lattice

interface use_walls
    module procedure :: use_walls_from_json
    module procedure :: use_walls_from_walls_potential
end interface use_walls

interface component_exists
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

interface components_interact
    module procedure :: components_interact_from_min_distance
    module procedure :: components_interact_from_pair_potential
    module procedure :: components_interact_from_ewald_real_pair
end interface components_interact

contains

    logical function apply_external_field_from_json(input_data, prefix) result(apply_external_field)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        apply_external_field = logical_from_json(input_data, prefix//"External Field.apply")
    end function apply_external_field_from_json

    logical function use_permittivity_from_json(input_data, prefix) result(use_permittivity)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_permittivity = logical_from_json(input_data, prefix//"Permittivity.use")
    end function

    pure logical function use_permittivity_from(permittivity) result(use_permittivity)
        class(Abstract_Permittivity), intent(in) :: permittivity

        select type (permittivity)
            type is (Concrete_Permittivity)
                use_permittivity = .true.
            class default
                use_permittivity = .false.
        end select
    end function

    logical function use_reciprocal_lattice_from_json(input_data, prefix) &
        result(use_reciprocal_lattice)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_reciprocal_lattice = logical_from_json(input_data, prefix//"Reciprocal Lattice.use")
    end function use_reciprocal_lattice_from_json

    pure logical function use_reciprocal_lattice_from(reciprocal_lattice) &
        result(use_reciprocal_lattice)
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice

        select type (reciprocal_lattice)
            type is (Concrete_Reciprocal_Lattice)
                use_reciprocal_lattice = .true.
            class default
                use_reciprocal_lattice = .false.
        end select
    end function use_reciprocal_lattice_from

    logical function use_walls_from_json(input_data, prefix) result(use_walls)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_walls = logical_from_json(input_data, prefix//"Walls.use")
    end function use_walls_from_json

    logical function use_walls_from_walls_potential(walls_potential) result(use_walls)
        class(Abstract_Walls_Potential), intent(in) :: walls_potential

        select type (walls_potential)
            type is (Concrete_Walls_Potential)
                use_walls = .true.
            class default
                use_walls = .false.
        end select
    end function use_walls_from_walls_potential

    pure logical function component_exists_from_number(component_number) result(component_exists)
        class(Abstract_Component_Number), intent(in) :: component_number

        select type (component_number)
            type is (Concrete_Component_Number)
                component_exists = .true.
            class default
                component_exists = .false.
        end select
    end function component_exists_from_number

    pure logical function components_interact_from_min_distance(min_distance) &
        result(components_interact)
        class(Abstract_Minimum_Distance), intent(in) :: min_distance

        select type (min_distance)
            type is (Concrete_Minimum_Distance)
                components_interact = .true.
            class default
                components_interact = .false.
        end select
    end function components_interact_from_min_distance

    pure logical function components_interact_from_pair_potential(pair_potential) &
        result(components_interact)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        select type (pair_potential)
            type is (Null_Pair_Potential)
                components_interact = .false.
            class default
                components_interact = .true.
        end select
    end function components_interact_from_pair_potential

    pure logical function components_interact_from_ewald_real_pair(ewald_real_pair) &
        result(components_interact)
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair

        select type (ewald_real_pair)
            type is (Null_Ewald_Real_Pair)
                components_interact = .false.
            class default
                components_interact = .true.
        end select
    end function components_interact_from_ewald_real_pair

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
        class(Abstract_Changed_Coordinates), intent(in) :: moved_positions

        select type (moved_positions)
            type is (Concrete_Moved_Positions)
                component_can_move = .true.
            class default
                component_can_move = .false.
        end select
    end function component_can_move

    pure logical function component_can_rotate(rotated_orientations)
        class(Abstract_Changed_Coordinates), intent(in) :: rotated_orientations

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                component_can_rotate = .true.
            class default
                component_can_rotate = .false.
        end select
    end function component_can_rotate

    logical function component_is_dipolar_from_json(input_data, prefix) &
        result(component_is_dipolar)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        component_is_dipolar = logical_from_json(input_data, prefix//"is dipolar")
    end function component_is_dipolar_from_json

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
            type is (Concrete_Component_Dipolar_Moments)
                component_is_dipolar = .true.
            class default
                component_is_dipolar = .false.
        end select
    end function component_is_dipolar_from_dipolar_moments

    logical function component_can_exchange_from_json(input_data, prefix) &
        result(component_can_exchange)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        component_can_exchange = logical_from_json(input_data, prefix//"can exchange")
    end function component_can_exchange_from_json

    pure logical function component_can_exchange_from_chemical_potential(&
        component_chemical_potential) result(component_can_exchange)
        class(Abstract_Component_Chemical_Potential), intent(in) :: component_chemical_potential

        select type (component_chemical_potential)
            type is (Concrete_Component_Chemical_Potential)
                component_can_exchange = .true.
            class default
                component_can_exchange = .false.
        end select
    end function component_can_exchange_from_chemical_potential

    pure logical function component_can_exchange_from_component_exchange(component_exchange) &
        result(component_can_exchange)
        class(Abstract_Component_Exchange), intent(in) :: component_exchange

        select type (component_exchange)
            type is (Concrete_Component_Exchange)
                component_can_exchange = .true.
            class default
                component_can_exchange = .false.
        end select
    end function component_can_exchange_from_component_exchange

    pure logical function component_can_change(moved_positions, rotated_orientations, &
        component_exchange)
        class(Abstract_Changed_Coordinates), intent(in) :: moved_positions, rotated_orientations
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