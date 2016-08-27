module procedures_moved_component_coordinates_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use procedures_mixture_inquirers, only: component_has_positions, component_has_orientations
use classes_moved_component_coordinates, only: Abstract_Moved_Component_Coordinates, &
    Null_Moved_Component_Coordinates
use classes_translated_positions, only: Concrete_Translated_Positions
use classes_rotated_orientations, only: Concrete_Rotated_Orientations
use module_move_tuning, only: Concrete_Move_Tuning_Parameters

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_translated_positions
    module procedure :: create_rotated_orientations
end interface create

contains

    subroutine create_translated_positions(translated_positions, periodic_box, positions, &
        tuning_parameters, generating_data, prefix)
        class(Abstract_Moved_Component_Coordinates), allocatable, intent(out) :: &
            translated_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: positions
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP), allocatable :: initial_delta(:)
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (component_has_positions(positions)) then
            allocate(Concrete_Translated_Positions :: translated_positions)
        else
            allocate(Null_Moved_Component_Coordinates :: translated_positions)
        end if

        select type (translated_positions)
            type is (Concrete_Translated_Positions)
                data_field = prefix//"Small Move.initial delta"
                call generating_data%get(data_field, initial_delta, data_found)
                call check_data_found(data_field, data_found)
                call translated_positions%construct(periodic_box, positions, initial_delta, &
                    tuning_parameters)
            type is (Null_Moved_Component_Coordinates)
                call translated_positions%construct()
            class default
                call error_exit("procedures_moved_component_coordinates_factory: "//&
                    "create_translated_positions: translated_positions: type unknown.")
        end select
    end subroutine create_translated_positions

    subroutine create_rotated_orientations(rotated_orientations, orientations, tuning_parameters, &
        generating_data, prefix)
        class(Abstract_Moved_Component_Coordinates), allocatable, intent(out) :: &
            rotated_orientations
        type(Concrete_Move_Tuning_Parameters), intent(in) :: tuning_parameters
        class(Abstract_Component_Coordinates), intent(in) :: orientations
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP) :: initial_delta
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (component_has_orientations(orientations)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Moved_Component_Coordinates :: rotated_orientations)
        end if

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                data_field = prefix//"Small Rotation.initial delta"
                call generating_data%get(data_field, initial_delta, data_found)
                call check_data_found(data_field, data_found)
                call rotated_orientations%construct(orientations, initial_delta, tuning_parameters)
            type is (Null_Moved_Component_Coordinates)
                call rotated_orientations%construct()
            class default
                call error_exit("create_rotated_orientations: rotated_orientations: type unknown")
        end select
    end subroutine create_rotated_orientations

    subroutine destroy(moved_coordinates)
        class(Abstract_Moved_Component_Coordinates), allocatable, intent(inout) :: moved_coordinates

        if (allocated(moved_coordinates)) then
            call moved_coordinates%destroy()
            deallocate(moved_coordinates)
        end if
    end subroutine destroy

end module procedures_moved_component_coordinates_factory
