module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity, Concrete_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice
use classes_external_field, only: Abstract_External_Field, Null_External_Field
use classes_visitable_walls, only: Abstract_Visitable_Walls, Concrete_Visitable_Walls
use classes_component_number, only: Abstract_Component_Number, Concrete_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations
use classes_min_distance, only: Abstract_Min_Distance, Concrete_Min_Distance
use classes_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Concrete_Component_Dipolar_Moments
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential
use classes_changed_coordinates, only: Abstract_Changed_Coordinates
use classes_moved_positions, only: Concrete_Moved_Positions
use classes_rotated_orientations, only: Concrete_Rotated_Orientations
use classes_component_exchange, only: Abstract_Component_Exchange, Concrete_Component_Exchange
use classes_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use classes_des_real_pair, only: Abstract_DES_Real_Pair, Null_DES_Real_Pair

implicit none

private
public :: periodicity_is_xyz, periodicity_is_xy, apply_external_field, use_permittivity, &
    use_reciprocal_lattice, use_walls, &
    component_exists, component_is_dipolar, component_has_positions, component_can_move, &
    component_has_orientations, component_can_exchange, component_can_rotate, components_interact, &
    component_interacts_with_wall, component_can_change

interface apply_external_field
    module procedure :: apply_external_field_from_json
    module procedure :: apply_external_field_from
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
    module procedure :: use_walls_from_walls
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
    module procedure :: components_interact_from_des_real_pair
end interface components_interact

contains

    pure logical function periodicity_is_xyz(periodic_box)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        select type (periodic_box)
            type is (XYZ_Periodic_Box)
                periodicity_is_xyz = .true.
            class default
                periodicity_is_xyz = .false.
        end select
    end function periodicity_is_xyz

    pure logical function periodicity_is_xy(periodic_box)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        select type (periodic_box)
            type is (XY_Periodic_Box)
                periodicity_is_xy = .true.
            class default
                periodicity_is_xy = .false.
        end select
    end function periodicity_is_xy

    logical function apply_external_field_from_json(input_data, prefix) result(apply_external_field)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        apply_external_field = logical_from_json(input_data, prefix//"External Field.apply")
    end function apply_external_field_from_json

    pure logical function apply_external_field_from(external_field) result(apply_external_field)
        class(Abstract_External_Field), intent(in) :: external_field

        select type (external_field)
            type is (Null_External_Field)
                apply_external_field = .false.
            class default
                apply_external_field = .true.
        end select
    end function apply_external_field_from

    logical function use_permittivity_from_json(input_data, prefix) result(use_permittivity)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        use_permittivity = logical_from_json(input_data, prefix//"Permittivity.use")
    end function use_permittivity_from_json

    pure logical function use_permittivity_from(permittivity) result(use_permittivity)
        class(Abstract_Permittivity), intent(in) :: permittivity

        select type (permittivity)
            type is (Concrete_Permittivity)
                use_permittivity = .true.
            class default
                use_permittivity = .false.
        end select
    end function use_permittivity_from

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

    logical function use_walls_from_walls(walls) result(use_walls)
        class(Abstract_Visitable_Walls), intent(in) :: walls

        select type (walls)
            type is (Concrete_Visitable_Walls)
                use_walls = .true.
            class default
                use_walls = .false.
        end select
    end function use_walls_from_walls

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
        class(Abstract_Min_Distance), intent(in) :: min_distance

        select type (min_distance)
            type is (Concrete_Min_Distance)
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

    pure logical function component_interacts_with_wall(wall_potential) result(component_interacts)
        class(Abstract_Pair_Potential), intent(in) :: wall_potential

        select type (wall_potential)
            type is (Null_Pair_Potential)
                component_interacts = .false.
            class default
                component_interacts = .true.
        end select
    end function component_interacts_with_wall

    pure logical function components_interact_from_des_real_pair(real_potential) &
        result(components_interact)
        class(Abstract_DES_Real_Pair), intent(in) :: real_potential

        select type (real_potential)
            type is (Null_DES_Real_Pair)
                components_interact = .false.
            class default
                components_interact = .true.
        end select
    end function components_interact_from_des_real_pair

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

    pure logical function component_is_dipolar_from_dipolar_moments(dipolar_moments) &
        result(component_is_dipolar)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments

        select type (dipolar_moments)
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

    logical function logical_from_json(input_data, statement)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: statement

        logical :: data_found

        call input_data%get(statement, logical_from_json, data_found)
        call check_data_found(statement, data_found)
    end function logical_from_json

end module procedures_property_inquirers
