module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box, XY_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity, Concrete_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice, Concrete_Reciprocal_Lattice
use classes_external_field, only: Abstract_External_Field, Null_External_Field
use classes_floor_penetration, only: Abstract_Floor_Penetration, Null_Floor_Penetration
use classes_visitable_walls, only: Abstract_Visitable_Walls, Concrete_Visitable_Walls
use classes_component_number, only: Abstract_Component_Number, Concrete_Component_Number
use classes_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations
use classes_min_distance, only: Abstract_Min_Distance, Concrete_Min_Distance
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments, &
    Concrete_Component_Dipole_Moments
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential
use classes_changed_box_size, only: Abstract_Changed_Box_Size, Null_Changed_Box_Size
use classes_moved_component_coordinates, only: Abstract_Moved_Component_Coordinates
use classes_translated_positions, only: Concrete_Translated_Positions
use classes_rotated_orientations, only: Concrete_Rotated_Orientations
use classes_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method, &
    Concrete_Particle_Insertion_Method
use classes_volume_change_method, only: Abstract_Volume_Change_Method, Concrete_Volume_Change_Method

implicit none

private
public :: periodicity_is_xyz, periodicity_is_xy, apply_external_field, use_permittivity, &
    use_reciprocal_lattice, use_walls, box_size_can_change, &
    component_exists, component_is_dipolar, component_has_positions, component_can_translate, &
    component_has_orientations, component_can_exchange, component_can_rotate, components_interact, &
    component_interacts_with_wall, measure_chemical_potentials, measure_pressure

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
    module procedure :: use_walls_from_floor_penetration
    module procedure :: use_walls_from_walls
end interface use_walls

interface box_size_can_change
    module procedure :: box_size_can_change_from_json
    module procedure :: box_size_can_change_from_changed_box_size
end interface box_size_can_change

interface component_exists
    module procedure :: component_exists_from_number
end interface component_exists

interface component_is_dipolar
    module procedure :: component_is_dipolar_from_json
    module procedure :: component_is_dipolar_from_dipole_moments
end interface component_is_dipolar

interface component_can_exchange
    module procedure :: component_can_exchange_from_json
    module procedure :: component_can_exchange_from_chemical_potential
end interface component_can_exchange

interface components_interact
    module procedure :: components_interact_from_min_distance
    module procedure :: components_interact_from_pair_potential
end interface components_interact

interface measure_chemical_potentials
    module procedure :: measure_chemical_potentials_from_json
    module procedure :: measure_chemical_potentials_from_particle_insertion_method
end interface measure_chemical_potentials

interface measure_pressure
    module procedure :: measure_pressure_from_json
    module procedure :: measure_pressure_from_volume_change_method
end interface measure_pressure

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

    logical function apply_external_field_from_json(generating_data, prefix) &
        result(apply_external_field)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        apply_external_field = logical_from_json(generating_data, prefix//"External Field.apply")
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

    logical function use_permittivity_from_json(generating_data, prefix) result(use_permittivity)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_permittivity = logical_from_json(generating_data, prefix//"Permittivity.use")
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

    logical function use_reciprocal_lattice_from_json(generating_data, prefix) &
        result(use_reciprocal_lattice)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_reciprocal_lattice = logical_from_json(generating_data, &
            prefix//"Reciprocal Lattice.use")
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

    logical function use_walls_from_json(generating_data, prefix) result(use_walls)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        use_walls = logical_from_json(generating_data, prefix//"Walls.use")
    end function use_walls_from_json

    pure logical function use_walls_from_floor_penetration(floor_penetration) result(use_walls)
        class(Abstract_Floor_Penetration), intent(in) :: floor_penetration

        select type (floor_penetration)
            type is (Null_Floor_Penetration)
                use_walls = .false.
            class default
                use_walls = .true.
        end select
    end function use_walls_from_floor_penetration

    pure logical function use_walls_from_walls(visitable_walls) result(use_walls)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        select type (visitable_walls)
            type is (Concrete_Visitable_Walls)
                use_walls = .true.
            class default
                use_walls = .false.
        end select
    end function use_walls_from_walls

    logical function box_size_can_change_from_json(generating_data, prefix) &
        result(box_size_can_change)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        box_size_can_change = logical_from_json(generating_data, prefix//"Box.size can change")
    end function box_size_can_change_from_json

    pure logical function box_size_can_change_from_changed_box_size(changed_box_size) &
        result(box_size_can_change)
        class(Abstract_Changed_Box_Size), intent(in) :: changed_box_size

        select type (changed_box_size)
            type is (Null_Changed_Box_Size)
                box_size_can_change = .false.
            class default
                box_size_can_change = .true.
        end select
    end function box_size_can_change_from_changed_box_size

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

    pure logical function component_has_positions(partcles_positions)
        class(Abstract_Component_Coordinates), intent(in) :: partcles_positions

        select type (partcles_positions)
            type is (Concrete_Component_Positions)
                component_has_positions = .true.
            class default
                component_has_positions = .false.
        end select
    end function component_has_positions

    pure logical function component_can_translate(translated_positions)
        class(Abstract_Moved_Component_Coordinates), intent(in) :: translated_positions

        select type (translated_positions)
            type is (Concrete_Translated_Positions)
                component_can_translate = .true.
            class default
                component_can_translate = .false.
        end select
    end function component_can_translate

    pure logical function component_can_rotate(rotated_orientations)
        class(Abstract_Moved_Component_Coordinates), intent(in) :: rotated_orientations

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                component_can_rotate = .true.
            class default
                component_can_rotate = .false.
        end select
    end function component_can_rotate

    logical function component_is_dipolar_from_json(generating_data, prefix) &
        result(component_is_dipolar)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        component_is_dipolar = logical_from_json(generating_data, prefix//"is dipolar")
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

    pure logical function component_is_dipolar_from_dipole_moments(dipole_moments) &
        result(component_is_dipolar)
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments

        select type (dipole_moments)
            type is (Concrete_Component_Dipole_Moments)
                component_is_dipolar = .true.
            class default
                component_is_dipolar = .false.
        end select
    end function component_is_dipolar_from_dipole_moments

    logical function component_can_exchange_from_json(generating_data, prefix) &
        result(component_can_exchange)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        component_can_exchange = logical_from_json(generating_data, prefix//"can exchange")
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

    logical function measure_chemical_potentials_from_json(exploring_data, prefix)&
        result(measure_chemical_potentials)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        measure_chemical_potentials = logical_from_json(exploring_data, &
            prefix//"mesure chemical potentials")
    end function measure_chemical_potentials_from_json

    pure logical function measure_chemical_potentials_from_particle_insertion_method&
        (particle_insertion_method) result(measure_chemical_potentials)
        class(Abstract_Particle_Insertion_Method), intent(in) :: particle_insertion_method

        select type (particle_insertion_method)
            type is (Concrete_Particle_Insertion_Method)
                measure_chemical_potentials = .true.
            class default
                measure_chemical_potentials = .false.
        end select
    end function measure_chemical_potentials_from_particle_insertion_method

    logical function measure_pressure_from_json(exploring_data, prefix) result(measure_pressure)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        measure_pressure = logical_from_json(exploring_data, prefix//"measure pressure")
    end function measure_pressure_from_json

    pure logical function measure_pressure_from_volume_change_method(volume_change_method) &
        result(measure_pressure)
        class(Abstract_Volume_Change_Method), intent(in) :: volume_change_method

        select type (volume_change_method)
            type is (Concrete_Volume_Change_Method)
                measure_pressure = .true.
            class default
                measure_pressure = .false.
        end select
    end function measure_pressure_from_volume_change_method

    logical function logical_from_json(input_data, statement)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: statement

        logical :: data_found

        call input_data%get(statement, logical_from_json, data_found)
        call check_data_found(statement, data_found)
    end function logical_from_json

end module procedures_property_inquirers
