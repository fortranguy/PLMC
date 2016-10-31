module procedures_mixture_inquirers

use json_module, only: json_file
use procedures_property_inquirers, only: logical_from_json
use procedures_checks, only: check_data_found, check_array_size
use classes_num_particles, only: Abstract_Num_Particles, Concrete_Num_Particles
use classes_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments, &
    Concrete_Component_Dipole_Moments
use classes_component_chemical_potential, only: Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential
use classes_moved_coordinates, only: Abstract_Moved_Coordinates
use classes_translated_positions, only: Concrete_Translated_Positions
use classes_rotated_orientations, only: Concrete_Rotated_Orientations

implicit none

private
public :: component_exists, component_is_dipolar, component_has_positions, &
    component_has_orientations, component_can_translate, component_can_rotate, &
    mixture_can_exchange, component_can_exchange, num_components, i_component, ij_components

interface component_exists
    module procedure :: component_exists_from_num
end interface component_exists

interface component_is_dipolar
    module procedure :: component_is_dipolar_from_json
    module procedure :: component_is_dipolar_from_dipole_moments
end interface component_is_dipolar

interface mixture_can_exchange
    module procedure :: mixture_can_exchange_from_json
end interface mixture_can_exchange

interface component_can_exchange
    module procedure :: component_can_exchange_from_chemical_potential
end interface component_can_exchange

interface num_components
    module procedure :: num_components_from_json
end interface num_components

interface i_component
    module procedure :: i_component_from_json
end interface i_component

interface ij_components
    module procedure :: ij_components_from_json
end interface ij_components

contains

    pure logical function component_exists_from_num(num_particles) result(component_exists)
        class(Abstract_Num_Particles), intent(in) :: num_particles

        select type (num_particles)
            type is (Concrete_Num_Particles)
                component_exists = .true.
            class default
                component_exists = .false.
        end select
    end function component_exists_from_num

    pure logical function component_has_positions(partcles_positions)
        class(Abstract_Component_Coordinates), intent(in) :: partcles_positions

        select type (partcles_positions)
            type is (Concrete_Component_Positions)
                component_has_positions = .true.
            class default
                component_has_positions = .false.
        end select
    end function component_has_positions

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

    pure logical function component_can_translate(translated_positions)
        class(Abstract_Moved_Coordinates), intent(in) :: translated_positions

        select type (translated_positions)
            type is (Concrete_Translated_Positions)
                component_can_translate = .true.
            class default
                component_can_translate = .false.
        end select
    end function component_can_translate

    pure logical function component_can_rotate(rotated_orientations)
        class(Abstract_Moved_Coordinates), intent(in) :: rotated_orientations

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                component_can_rotate = .true.
            class default
                component_can_rotate = .false.
        end select
    end function component_can_rotate

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

    logical function mixture_can_exchange_from_json(generating_data, prefix) &
        result(mixture_can_exchange)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        mixture_can_exchange = logical_from_json(generating_data, &
            prefix//"can exchange with reservoir")
    end function mixture_can_exchange_from_json

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

    integer function num_components_from_json(generating_data, prefix) result(num_components)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"number of components"
        call generating_data%get(data_field, num_components, data_found)
        call check_data_found(data_field, data_found)
    end function num_components_from_json

    integer function i_component_from_json(exploring_data, prefix) result(i_component)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"component number"
        call exploring_data%get(data_field, i_component, data_found)
        call check_data_found(data_field, data_found)
    end function i_component_from_json

    function ij_components_from_json(exploring_data, prefix) result(ij_components)
        integer :: ij_components(2)
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        integer, allocatable :: json_ij_components(:)
        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//"components couple"
        call exploring_data%get(data_field, json_ij_components, data_found)
        call check_data_found(data_field, data_found)
        call check_array_size("procedures_property_inquirers: ij_components_from_json", &
            "json_ij_components", json_ij_components, size(ij_components))
        ij_components = json_ij_components
    end function ij_components_from_json

end module procedures_mixture_inquirers
