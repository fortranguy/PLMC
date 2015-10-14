module procedures_property_inquirers

use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_walls_potential, only: Abstract_Walls_Potential, Null_Walls_Potential
use class_particles_number, only: Abstract_Particles_Number, Null_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter, Null_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, Null_Particles_Moment_Norm
use class_component_positions, only: Abstract_Component_Positions, Null_Component_Positions
use class_component_orientations, only: Abstract_Component_Orientations, Null_Component_Orientations
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments, &
    Null_Particles_Dipolar_Moments
use class_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Null_Component_Chemical_Potential
use class_moved_positions, only: Abstract_Moved_Positions, Null_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations, Null_Rotated_Orientations
use class_particles_exchange, only: Abstract_Particles_Exchange, Null_Particles_Exchange
use class_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, Null_Ewald_Real_Pair

implicit none

private
public :: apply_external_field, use_reciprocal_lattice, use_walls, &
    particles_exist, particles_have_positions, particles_can_move, particles_are_dipolar, &
    particles_have_orientations, particles_can_exchange, particles_can_rotate, particles_interact, &
    particles_can_change

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

interface particles_exist
    module procedure :: particles_exist_from_json
    module procedure :: particles_exist_from_number
end interface particles_exist

interface particles_are_dipolar
    module procedure :: particles_are_dipolar_from_json
    module procedure :: particles_are_dipolar_from_dipolar_moments
end interface particles_are_dipolar

interface particles_can_exchange
    module procedure :: particles_can_exchange_from_json
    module procedure :: particles_can_exchange_from_chemical_potential
    module procedure :: particles_can_exchange_from_particles_exchange
end interface particles_can_exchange

interface particles_interact
    module procedure :: particles_interact_short
    module procedure :: particles_interact_long
end interface particles_interact

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

    logical function particles_exist_from_json(input_data, prefix) result(particles_exist)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        particles_exist = logical_from_json(input_data, prefix//"exist")
    end function particles_exist_from_json

    pure logical function particles_exist_from_number(particles_number) result(particles_exist)
        class(Abstract_Particles_Number), intent(in) :: particles_number

        select type (particles_number)
            type is (Null_Particles_Number)
                particles_exist = .false.
            class default
                particles_exist = .true.
        end select
    end function particles_exist_from_number

    pure logical function particles_interact_short(pair_potential) result(particles_interact)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        select type (pair_potential)
            type is (Null_Pair_Potential)
                particles_interact = .false.
            class default
                particles_interact = .true.
        end select
    end function particles_interact_short

    pure logical function particles_interact_long(ewald_pair) result(particles_interact)
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_pair

        select type (ewald_pair)
            type is (Null_Ewald_Real_Pair)
                particles_interact = .false.
            class default
                particles_interact = .true.
        end select
    end function particles_interact_long

    pure logical function particles_have_positions(partcles_positions)
        class(Abstract_Component_Positions), intent(in) :: partcles_positions

        select type (partcles_positions)
            type is (Null_Component_Positions)
                particles_have_positions = .false.
            class default
                particles_have_positions = .true.
        end select
    end function particles_have_positions

    pure logical function particles_can_move(moved_positions)
        class(Abstract_Moved_Positions), intent(in) :: moved_positions

        select type (moved_positions)
            type is (Null_Moved_Positions)
                particles_can_move = .false.
            class default
                particles_can_move = .true.
        end select
    end function particles_can_move

    pure logical function particles_can_rotate(rotated_orientations)
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations

        select type (rotated_orientations)
            type is (Null_Rotated_Orientations)
                particles_can_rotate = .false.
            class default
                particles_can_rotate = .true.
        end select
    end function particles_can_rotate

    logical function particles_are_dipolar_from_json(input_data, prefix) &
        result(particles_are_dipolar)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (particles_exist(input_data, prefix)) then
            particles_are_dipolar = logical_from_json(input_data, prefix//"are dipolar")
        else
            particles_are_dipolar = .false.
        end if
    end function particles_are_dipolar_from_json

    pure logical function particles_are_dipolar_from_moment_norm(particles_moment_norm) &
        result(particles_are_dipolar)
        class(Abstract_Particles_Moment_Norm), intent(in) :: particles_moment_norm

        select type (particles_moment_norm)
            type is (Null_Particles_Moment_Norm)
                particles_are_dipolar = .false.
            class default
                particles_are_dipolar = .true.
        end select
    end function particles_are_dipolar_from_moment_norm

    pure logical function particles_have_orientations(particles_orientations)
        class(Abstract_Component_Orientations), intent(in) :: particles_orientations

        select type (particles_orientations)
            type is (Null_Component_Orientations)
                particles_have_orientations = .false.
            class default
                particles_have_orientations = .true.
        end select
    end function particles_have_orientations

    pure logical function particles_are_dipolar_from_dipolar_moments(particles_dipolar_moments) &
        result(particles_are_dipolar)
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        select type (particles_dipolar_moments)
            type is (Null_Particles_Dipolar_Moments)
                particles_are_dipolar = .false.
            class default
                particles_are_dipolar = .true.
        end select
    end function particles_are_dipolar_from_dipolar_moments

    logical function particles_can_exchange_from_json(input_data, prefix) &
        result(particles_can_exchange)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        if (particles_exist(input_data, prefix)) then
            particles_can_exchange = logical_from_json(input_data, prefix//"can exchange")
        else
            particles_can_exchange = .false.
        end if
    end function particles_can_exchange_from_json

    pure logical function particles_can_exchange_from_chemical_potential(&
        particles_chemical_potential) result(particles_can_exchange)
        class(Abstract_Component_Chemical_Potential), intent(in) :: particles_chemical_potential

        select type (particles_chemical_potential)
            type is (Null_Component_Chemical_Potential)
                particles_can_exchange = .false.
            class default
                particles_can_exchange = .true.
        end select
    end function particles_can_exchange_from_chemical_potential

    pure logical function particles_can_exchange_from_particles_exchange(particles_exchange) &
        result(particles_can_exchange)
        class(Abstract_Particles_Exchange), intent(in) :: particles_exchange

        select type (particles_exchange)
            type is (Null_Particles_Exchange)
                particles_can_exchange = .false.
            class default
                particles_can_exchange = .true.
        end select
    end function particles_can_exchange_from_particles_exchange

    pure logical function particles_can_change(moved_positions, rotated_orientations, &
        particles_exchange)
        class(Abstract_Moved_Positions), intent(in) :: moved_positions
        class(Abstract_Rotated_Orientations), intent(in) :: rotated_orientations
        class(Abstract_Particles_Exchange), intent(in) :: particles_exchange

        particles_can_change = particles_can_move(moved_positions) .or. &
            particles_can_rotate(rotated_orientations) .or. &
            particles_can_exchange(particles_exchange)
    end function particles_can_change

    logical function logical_from_json(input_data, switch_name)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: switch_name

        logical :: data_found

        call input_data%get(switch_name, logical_from_json, data_found)
        call check_data_found(switch_name, data_found)
    end function logical_from_json

end module procedures_property_inquirers
