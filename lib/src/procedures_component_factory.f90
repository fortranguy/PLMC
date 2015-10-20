module procedures_component_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use procedures_coordinates_micro, only: read_coordinates
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_number, only: Abstract_Component_Number, Concrete_Component_Number, &
    Null_Component_Number
use class_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations, Null_Component_Coordinates
use class_component_chemical_potential, only : Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential, Null_Component_Chemical_Potential
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Concrete_Component_Dipolar_Moments, Null_Component_Dipolar_Moments
use class_component_total_moment, only: Abstract_Component_Total_Moment, &
    Concrete_Component_Total_Moment, Null_Component_Total_Moment
use types_component_wrapper, only: Component_Wrapper
use procedures_property_inquirers, only: use_walls, component_is_dipolar, &
    component_can_exchange

implicit none

private
public :: component_create, component_set, component_destroy

interface component_create
    module procedure :: create_all
    module procedure :: create_number
    module procedure :: create_positions
    module procedure :: create_orientations
    module procedure :: create_dipolar_moments
    module procedure :: create_total_moment
    module procedure :: create_chemical_potential
end interface component_create

interface component_set
    module procedure :: set_coordinates
end interface component_set

interface component_destroy
    module procedure :: destroy_chemical_potential
    module procedure :: destroy_total_moment
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

        call component_create(component%number, exists, input_data, prefix)
        call component_create(component%positions, exists, periodic_box, component%number)
        call component_set(component%positions, input_data, prefix//"initial positions")
        is_dipolar = exists .and. component_is_dipolar(input_data, prefix)
        call component_create(component%orientations, is_dipolar, component%number)
        call component_set(component%orientations, input_data, &
            prefix//"initial orientations")
        call component_create(component%dipolar_moments, is_dipolar, component%orientations, &
            input_data, prefix)
        call component_create(component%total_moment, component%dipolar_moments)
        can_exchange = exists .and. component_can_exchange(input_data, prefix)
        call component_create(component%chemical_potential, can_exchange, input_data, &
            prefix)
    end subroutine create_all

    subroutine destroy_all(component)
        type(Component_Wrapper), intent(inout) :: component

        call component_destroy(component%chemical_potential)
        call component_destroy(component%total_moment)
        call component_destroy(component%dipolar_moments)
        call component_destroy(component%orientations)
        call component_destroy(component%positions)
        call component_destroy(component%number)
    end subroutine destroy_all

    subroutine create_number(component_number, exists, input_data, prefix)
        class(Abstract_Component_Number), allocatable, intent(out) :: component_number
        logical, intent(in) :: exists
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_particles

        if (exists) then
            allocate(Concrete_Component_Number :: component_number)
            data_field = prefix//"number"
            call input_data%get(data_field, num_particles, data_found)
            call check_data_found(data_field, data_found)
            deallocate(data_field)
        else
            num_particles = 0
            allocate(Null_Component_Number :: component_number)
        end if
        call component_number%set(num_particles)
    end subroutine create_number

    subroutine destroy_number(component_number)
        class(Abstract_Component_Number), allocatable, intent(inout) :: component_number

        if (allocated(component_number)) deallocate(component_number)
    end subroutine destroy_number

    subroutine create_positions(component_positions, exists, periodic_box, component_number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: component_positions
        logical, intent(in) :: exists
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Number), intent(in) :: component_number

        if (exists) then
            allocate(Concrete_Component_Positions :: component_positions)
        else
            allocate(Null_Component_Coordinates :: component_positions)
        end if
        select type (component_positions)
            type is (Concrete_Component_Positions)
                call component_positions%construct(periodic_box, component_number)
            type is (Null_Component_Coordinates)
                call component_positions%construct()
        end select
    end subroutine create_positions

    subroutine create_orientations(component_orientations, is_dipolar, component_number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: component_orientations
        logical, intent(in) :: is_dipolar
        class(Abstract_Component_Number), intent(in) :: component_number

        if (is_dipolar) then
            allocate(Concrete_Component_Orientations :: component_orientations)
        else
            allocate(Null_Component_Coordinates :: component_orientations)
        end if
        select type (component_orientations)
            type is (Concrete_Component_Orientations)
                call component_orientations%construct(component_number)
            type is (Null_Component_Coordinates)
                call component_orientations%destroy()
        end select
    end subroutine create_orientations

    subroutine destroy_coordinates(component_coordinates)
        class(Abstract_Component_Coordinates), allocatable, intent(inout) :: component_coordinates

        call component_coordinates%destroy()
        if (allocated(component_coordinates)) deallocate(component_coordinates)
    end subroutine destroy_coordinates

    subroutine set_coordinates(component_coordinates, input_data, data_field)
        class(Abstract_Component_Coordinates), intent(inout) :: component_coordinates
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: data_field

        character(len=:), allocatable :: filename
        logical :: data_found
        real(DP), allocatable :: file_coordinates(:, :)
        integer :: i_particle

        if (component_coordinates%get_num() == 0) return
        call input_data%get(data_field, filename, data_found)
        call check_data_found(data_field, data_found)
        call read_coordinates(file_coordinates, filename)
        if (size(file_coordinates, 2) /= component_coordinates%get_num()) then
            call error_exit("set_coordinates from "//filename//": wrong number of lines.")
        end if
        do i_particle = 1, component_coordinates%get_num()
            call component_coordinates%set(i_particle, file_coordinates(:, i_particle))
        end do
        if (allocated(file_coordinates)) deallocate(file_coordinates)
        deallocate(filename)
    end subroutine set_coordinates

    subroutine create_dipolar_moments(component_dipolar_moments, is_dipolar, &
        component_orientations, input_data, prefix)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(out) :: &
            component_dipolar_moments
        logical, intent(in) :: is_dipolar
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        if (is_dipolar) then
            data_field = prefix//"moment norm"
            call input_data%get(data_field, moment_norm, data_found)
            call check_data_found(data_field, data_found)
            deallocate(data_field)
            allocate(Concrete_Component_Dipolar_Moments :: component_dipolar_moments)
        else
            moment_norm = 0._DP
            allocate(Null_Component_Dipolar_Moments :: component_dipolar_moments)
        end if
        call component_dipolar_moments%construct(moment_norm, component_orientations)
    end subroutine create_dipolar_moments

    subroutine destroy_dipolar_moments(component_dipolar_moments)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(inout) :: &
            component_dipolar_moments

        call component_dipolar_moments%destroy()
        if (allocated(component_dipolar_moments)) deallocate(component_dipolar_moments)
    end subroutine destroy_dipolar_moments

    subroutine create_total_moment(component_total_moment, component_dipolar_moments)
        class(Abstract_Component_Total_Moment), allocatable, intent(out) :: component_total_moment
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments

        if (component_is_dipolar(component_dipolar_moments)) then
            allocate(Concrete_Component_Total_Moment :: component_total_moment)
        else
            allocate(Null_Component_Total_Moment :: component_total_moment)
        end if
        call component_total_moment%construct(component_dipolar_moments)
    end subroutine create_total_moment

    subroutine destroy_total_moment(component_total_moment)
        class(Abstract_Component_Total_Moment), allocatable, intent(inout) :: component_total_moment

        call component_total_moment%destroy()
        if (allocated(component_total_moment)) deallocate(component_total_moment)
    end subroutine destroy_total_moment

    subroutine create_chemical_potential(component_chemical_potential, can_exchange, input_data, &
            prefix)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(out) :: &
            component_chemical_potential
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
            allocate(Concrete_Component_Chemical_Potential :: component_chemical_potential)
            deallocate(data_field)
        else
            allocate(Null_Component_Chemical_Potential :: component_chemical_potential)
        end if
        call component_chemical_potential%set(density, excess)
    end subroutine create_chemical_potential

    subroutine destroy_chemical_potential(component_chemical_potential)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(inout) :: &
            component_chemical_potential

        if (allocated(component_chemical_potential)) deallocate(component_chemical_potential)
    end subroutine destroy_chemical_potential

end module procedures_component_factory
