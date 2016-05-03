module procedures_changed_coordinates_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_changed_coordinates, only: Abstract_Changed_Coordinates, Null_Changed_Coordinates
use classes_moved_positions, only: Concrete_Moved_Positions
use classes_rotated_orientations, only: Concrete_Rotated_Orientations
use module_change_tuning, only: Concrete_Change_Tuning_Parameters
use procedures_property_inquirers, only: component_has_positions, component_has_orientations

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_moved_positions
    module procedure :: create_rotated_orientations
end interface create

contains

    subroutine create_moved_positions(moved_positions, periodic_box, positions, tuning_parameters, &
        generating_data, prefix)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: moved_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: positions
        type(Concrete_Change_Tuning_Parameters), intent(in) :: tuning_parameters
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP), allocatable :: delta(:)
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (component_has_positions(positions)) then
            allocate(Concrete_Moved_Positions :: moved_positions)
        else
            allocate(Null_Changed_Coordinates :: moved_positions)
        end if

        select type (moved_positions)
            type is (Concrete_Moved_Positions)
                data_field = prefix//"Small Move.initial delta"
                call generating_data%get(data_field, delta, data_found)
                call check_data_found(data_field, data_found)
                call moved_positions%construct(periodic_box, positions, delta, tuning_parameters)
            type is (Null_Changed_Coordinates)
                call moved_positions%construct()
            class default
                call error_exit("create_moved_positions: moved_positions: type unknown")
        end select
    end subroutine create_moved_positions

    subroutine create_rotated_orientations(rotated_orientations, orientations, tuning_parameters, &
        generating_data, prefix)
        class(Abstract_Changed_Coordinates), allocatable, intent(out) :: rotated_orientations
        type(Concrete_Change_Tuning_Parameters), intent(in) :: tuning_parameters
        class(Abstract_Component_Coordinates), intent(in) :: orientations
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        real(DP) :: delta
        character(len=:), allocatable :: data_field
        logical :: data_found

        if (component_has_orientations(orientations)) then
            allocate(Concrete_Rotated_Orientations :: rotated_orientations)
        else
            allocate(Null_Changed_Coordinates :: rotated_orientations)
        end if

        select type (rotated_orientations)
            type is (Concrete_Rotated_Orientations)
                data_field = prefix//"Small Rotation.initial delta"
                call generating_data%get(data_field, delta, data_found)
                call check_data_found(data_field, data_found)
                call rotated_orientations%construct(orientations, delta, tuning_parameters)
            type is (Null_Changed_Coordinates)
                call rotated_orientations%construct()
            class default
                call error_exit("create_rotated_orientations: rotated_orientations: type unknown")
        end select
    end subroutine create_rotated_orientations

    subroutine destroy(changed_coordinates)
        class(Abstract_Changed_Coordinates), allocatable, intent(inout) :: changed_coordinates

        if (allocated(changed_coordinates)) then
            call changed_coordinates%destroy()
            deallocate(changed_coordinates)
        end if
    end subroutine destroy

end module procedures_changed_coordinates_factory
