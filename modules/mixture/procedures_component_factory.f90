module procedures_component_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_number, only: Abstract_Component_Number, Concrete_Component_Number, &
    Null_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations, Null_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Concrete_Component_Dipolar_Moments, Null_Component_Dipolar_Moments
use class_component_chemical_potential, only : Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential, Null_Component_Chemical_Potential
use class_component_average_number, only: Abstract_Component_Average_Number, &
    Constant_Component_Average_Number, Variable_Component_Average_Number, &
    Null_Component_Average_Number
use types_component_wrapper, only: Component_Wrapper
use procedures_property_inquirers, only: use_walls, component_is_dipolar, component_can_exchange

implicit none

private
public :: component_create, component_destroy

interface component_create
    module procedure :: create_all
    module procedure :: create_number
    module procedure :: create_positions
    module procedure :: create_orientations
    module procedure :: create_dipolar_moments
    module procedure :: create_chemical_potential
    module procedure :: create_average_number
end interface component_create

interface component_destroy
    module procedure :: destroy_average_number
    module procedure :: destroy_chemical_potential
    module procedure :: destroy_dipolar_moments
    module procedure :: destroy_coordinates
    module procedure :: destroy_number
    module procedure :: destroy_all
end interface component_destroy

contains

    subroutine create_all(component, exists, periodic_box, input_data, prefix)
        type(Component_Wrapper), intent(out) :: component
        logical, intent(in) :: exists
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: is_dipolar, can_exchange

        call component_create(component%number, exists)
        call component_create(component%positions, exists, periodic_box, component%number)
        is_dipolar = exists .and. component_is_dipolar(input_data, prefix)
        call component_create(component%orientations, is_dipolar, component%number)
        call component_create(component%dipolar_moments, is_dipolar, component%orientations, &
            input_data, prefix)
        can_exchange = exists .and. component_can_exchange(input_data, prefix)
        call component_create(component%chemical_potential, can_exchange, input_data, prefix)
        call component_create(component%average_number, periodic_box, component%number, component%&
            chemical_potential)
    end subroutine create_all

    subroutine destroy_all(component)
        type(Component_Wrapper), intent(inout) :: component

        call component_destroy(component%average_number)
        call component_destroy(component%chemical_potential)
        call component_destroy(component%dipolar_moments)
        call component_destroy(component%orientations)
        call component_destroy(component%positions)
        call component_destroy(component%number)
    end subroutine destroy_all

    !> Number will be set with coordinates, cf. [[Abstract_Coordinates_Reader]]. Too fragile?
    subroutine create_number(number, exists)
        class(Abstract_Component_Number), allocatable, intent(out) :: number
        logical, intent(in) :: exists

        if (exists) then
            allocate(Concrete_Component_Number :: number)
        else
            allocate(Null_Component_Number :: number)
        end if
    end subroutine create_number

    subroutine destroy_number(number)
        class(Abstract_Component_Number), allocatable, intent(inout) :: number

        if (allocated(number)) deallocate(number)
    end subroutine destroy_number

    subroutine create_positions(positions, exists, periodic_box, number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: positions
        logical, intent(in) :: exists
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Number), intent(in) :: number

        if (exists) then
            allocate(Concrete_Component_Positions :: positions)
        else
            allocate(Null_Component_Coordinates :: positions)
        end if
        select type (positions)
            type is (Concrete_Component_Positions)
                call positions%construct(periodic_box, number)
            type is (Null_Component_Coordinates)
                call positions%construct()
            class default
                call error_exit("create_positions: positions: unknown type.")
        end select
    end subroutine create_positions

    subroutine create_orientations(orientations, is_dipolar, number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: orientations
        logical, intent(in) :: is_dipolar
        class(Abstract_Component_Number), intent(in) :: number

        if (is_dipolar) then
            allocate(Concrete_Component_Orientations :: orientations)
        else
            allocate(Null_Component_Coordinates :: orientations)
        end if
        select type (orientations)
            type is (Concrete_Component_Orientations)
                call orientations%construct(number)
            type is (Null_Component_Coordinates)
                call orientations%destroy()
            class default
                call error_exit("create_orientations: orientations: unknown type.")
        end select
    end subroutine create_orientations

    subroutine destroy_coordinates(coordinates)
        class(Abstract_Component_Coordinates), allocatable, intent(inout) :: coordinates

        if (allocated(coordinates)) then
            call coordinates%destroy()
            deallocate(coordinates)
        end if
    end subroutine destroy_coordinates

    subroutine create_dipolar_moments(dipolar_moments, is_dipolar, orientations, input_data, prefix)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(out) :: dipolar_moments
        logical, intent(in) :: is_dipolar
        class(Abstract_Component_Coordinates), intent(in) :: orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        if (is_dipolar) then
            data_field = prefix//"moment norm"
            call input_data%get(data_field, moment_norm, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Dipolar_Moments :: dipolar_moments)
        else
            moment_norm = 0._DP
            allocate(Null_Component_Dipolar_Moments :: dipolar_moments)
        end if
        call dipolar_moments%construct(moment_norm, orientations)
    end subroutine create_dipolar_moments

    subroutine destroy_dipolar_moments(dipolar_moments)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(inout) :: dipolar_moments

        if (allocated(dipolar_moments)) then
            call dipolar_moments%destroy()
            deallocate(dipolar_moments)
        end if
    end subroutine destroy_dipolar_moments

    subroutine create_chemical_potential(chemical_potential, can_exchange, input_data, prefix)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(out) :: &
            chemical_potential
        logical, intent(in) :: can_exchange
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, excess

        if (can_exchange) then
            data_field = prefix//"Chemical Potential.density"
            call input_data%get(data_field, density, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"Chemical Potential.excess"
            call input_data%get(data_field, excess, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Chemical_Potential :: chemical_potential)
        else
            allocate(Null_Component_Chemical_Potential :: chemical_potential)
        end if
        call chemical_potential%set(density, excess)
    end subroutine create_chemical_potential

    subroutine destroy_chemical_potential(chemical_potential)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(inout) :: &
            chemical_potential

        if (allocated(chemical_potential)) deallocate(chemical_potential)
    end subroutine destroy_chemical_potential

    subroutine create_average_number(average_number, periodic_box, number, chemical_potential)
        class(Abstract_Component_Average_Number), allocatable, intent(out) :: average_number
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Number), intent(in) :: number
        class(Abstract_Component_Chemical_Potential), intent(in) :: chemical_potential

        if (component_can_exchange(chemical_potential)) then
            allocate(Variable_Component_Average_Number :: average_number)
        else
            allocate(Constant_Component_Average_Number :: average_number)
        end if

        select type (average_number)
            type is (Constant_Component_Average_Number)
                call average_number%construct(number)
            type is (Variable_Component_Average_Number)
                call average_number%construct(periodic_box, chemical_potential)
            type is (Null_Component_Average_Number)
                call average_number%construct()
            class default
                call error_exit("create_average_number: average_number: type unknown.")
        end select
    end subroutine create_average_number

    subroutine destroy_average_number(average_number)
        class(Abstract_Component_Average_Number), allocatable, intent(inout) :: average_number

        if (allocated(average_number)) then
            call average_number%destroy()
            deallocate(average_number)
        end if
    end subroutine destroy_average_number

end module procedures_component_factory
