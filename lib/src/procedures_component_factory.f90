module procedures_component_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_errors, only: error_exit
use procedures_checks, only: check_data_found
use procedures_coordinates_micro, only: read_coordinates
use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use class_component_number, only: Abstract_Component_Number, &
    Concrete_Component_Number, Null_Component_Number
use class_component_diameter, only: Abstract_Component_Diameter, &
    Concrete_Component_Diameter, Null_Component_Diameter
use class_component_moment_norm, only: Abstract_Component_Moment_Norm, &
    Concrete_Component_Moment_Norm, Null_Component_Moment_Norm
use class_component_coordinates, only: Abstract_Component_Coordinates, &
    Concrete_Component_Positions, Concrete_Component_Orientations, Null_Component_Coordinates
use class_component_chemical_potential, only : Abstract_Component_Chemical_Potential, &
    Concrete_Component_Chemical_Potential, Null_Component_Chemical_Potential
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments, &
    Concrete_Component_Dipolar_Moments, Null_Component_Dipolar_Moments
use class_component_total_moment, only: Abstract_Component_Total_Moment, &
    Concrete_Component_Total_Moment, Null_Component_Total_Moment
use types_component_wrapper, only: Component_Wrapper, Mixture_Wrapper_Old
use procedures_property_inquirers, only: use_walls, component_exists, component_is_dipolar, &
    component_can_exchange

implicit none

private
public :: component_factory_create, component_factory_set, component_factory_destroy

interface component_factory_create
    module procedure :: component_factory_create_all
    module procedure :: allocate_and_set_number
    module procedure :: allocate_and_set_diameter
    module procedure :: allocate_and_set_inter_diameter
    module procedure :: allocate_and_set_moment_norm
    module procedure :: allocate_and_construct_positions
    module procedure :: allocate_and_construct_orientations
    module procedure :: allocate_and_construct_dipolar_moments
    module procedure :: allocate_and_construct_total_moment
    module procedure :: allocate_and_set_chemical_potential
end interface component_factory_create

interface component_factory_set
    module procedure :: set_coordinates
end interface component_factory_set

interface component_factory_destroy
    module procedure :: deallocate_chemical_potential
    module procedure :: destroy_and_deallocate_total_moment
    module procedure :: destroy_and_deallocate_dipolar_moments
    module procedure :: destroy_and_deallocate_coordinates
    module procedure :: deallocate_moment_norm
    module procedure :: deallocate_diameter
    module procedure :: deallocate_number
    module procedure :: component_factory_destroy_all
end interface component_factory_destroy

contains

    subroutine component_factory_create_all(component, environment, input_data, prefix)
        type(Component_Wrapper), intent(out) :: component
        type(Environment_Wrapper), intent(in) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: exist, exists_and_walls_used, dipolar

        exist = component_exists(input_data, prefix)
        call component_factory_create(component%number, exist, input_data, prefix)
        call component_factory_create(component%diameter, exist, input_data, prefix)
        exists_and_walls_used = exist .and. use_walls(environment%walls_potential)
        call component_factory_create(component%wall_diameter, exists_and_walls_used, input_data, &
            prefix//"With Walls.")
        dipolar = component_is_dipolar(input_data, prefix)
        call component_factory_create(component%moment_norm, dipolar, input_data, prefix)
        call component_factory_create(component%positions, exist, environment%periodic_box, &
            component%number)
        call component_factory_set(component%positions, input_data, prefix//"initial positions")
        call component_factory_create(component%orientations, dipolar, component%number)
        call component_factory_set(component%orientations, input_data, prefix//"initial orientations")
        call component_factory_create(component%dipolar_moments, dipolar, component%moment_norm, &
            component%orientations)
        call component_factory_create(component%total_moment, component%dipolar_moments)
        call component_factory_create(component%chemical_potential, input_data, prefix)
    end subroutine component_factory_create_all

    subroutine allocate_and_set_number(component_number, exist, input_data, prefix)
        class(Abstract_Component_Number), allocatable, intent(out) :: component_number
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_component ! ambiguous

        if (exist) then
            data_field = prefix//"number"
            call input_data%get(data_field, num_component, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Number :: component_number)
            deallocate(data_field)
        else
            allocate(Null_Component_Number :: component_number)
        end if
        call component_number%set(num_component)
    end subroutine allocate_and_set_number

    subroutine allocate_and_set_diameter(component_diameter, exist, input_data, prefix)
        class(Abstract_Component_Diameter), allocatable, intent(out) :: component_diameter
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_diameter(component_diameter, exist)
        call set_diameter(component_diameter, exist, input_data, prefix)
    end subroutine allocate_and_set_diameter

    subroutine allocate_diameter(component_diameter, exist)
        class(Abstract_Component_Diameter), allocatable, intent(out) :: component_diameter
        logical, intent(in) :: exist

        if (exist) then
            allocate(Concrete_Component_Diameter :: component_diameter)
        else
            allocate(Null_Component_Diameter :: component_diameter)
        end if
    end subroutine allocate_diameter

    subroutine set_diameter(component_diameter, exist, input_data, prefix)
        class(Abstract_Component_Diameter), intent(inout) :: component_diameter
        logical, intent(in) :: exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: diameter, diameter_min_factor

        if (exist) then
            data_field = prefix//"diameter"
            call input_data%get(data_field, diameter, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"minimum diameter factor"
            call input_data%get(data_field, diameter_min_factor, data_found)
            call check_data_found(data_field, data_found)
            deallocate(data_field)
        else
            diameter = 0._DP
            diameter_min_factor = 0._DP
        end if
        call component_diameter%set(diameter, diameter_min_factor)
    end subroutine set_diameter

    subroutine allocate_and_set_inter_diameter(inter_component_diameter, exists, &
        component_diameter_1, component_diameter_2, input_data, prefix)
        class(Abstract_Component_Diameter), allocatable, intent(out) :: inter_component_diameter
        logical, intent(in) :: exists
        class(Abstract_Component_Diameter), intent(in) :: component_diameter_1, component_diameter_2
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: inter_diameter, inter_diameter_min_factor
        real(DP) :: inter_diameter_offset

        if (exists) then
            data_field = prefix//"offset"
            call input_data%get(data_field, inter_diameter_offset, data_found)
            call check_data_found(data_field, data_found)
            data_field = prefix//"minimum diameter factor"
            call input_data%get(data_field, inter_diameter_min_factor, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Diameter :: inter_component_diameter)
            deallocate(data_field)
        else
            inter_diameter_offset = 0._DP
            allocate(Null_Component_Diameter :: inter_component_diameter)
        end if
        inter_diameter = (component_diameter_1%get() + component_diameter_2%get()) / 2._DP + &
            inter_diameter_offset
        call inter_component_diameter%set(inter_diameter, inter_diameter_min_factor)
    end subroutine allocate_and_set_inter_diameter

    subroutine allocate_and_set_moment_norm(component_moment_norm, dipolar, input_data, prefix)
        class(Abstract_Component_Moment_Norm), allocatable, intent(out) :: component_moment_norm
        logical, intent(in) :: dipolar
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        if (dipolar) then
            data_field = prefix//"moment norm"
            call input_data%get(data_field, moment_norm, data_found)
            call check_data_found(data_field, data_found)
            allocate(Concrete_Component_Moment_Norm :: component_moment_norm)
            deallocate(data_field)
        else
            moment_norm = 0._DP
            allocate(Null_Component_Moment_Norm :: component_moment_norm)
        end if
        call component_moment_norm%set(moment_norm)
    end subroutine allocate_and_set_moment_norm

    subroutine allocate_and_construct_positions(component_positions, exist, periodic_box, &
        component_number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: component_positions
        logical, intent(in) :: exist
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Number), intent(in) :: component_number

        if (exist) then
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
    end subroutine allocate_and_construct_positions

    subroutine allocate_and_construct_orientations(component_orientations, dipolar, &
        component_number)
        class(Abstract_Component_Coordinates), allocatable, intent(out) :: component_orientations
        logical, intent(in) :: dipolar
        class(Abstract_Component_Number), intent(in) :: component_number

        if (dipolar) then
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
    end subroutine allocate_and_construct_orientations

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

    subroutine allocate_and_construct_dipolar_moments(component_dipolar_moments, dipolar, &
        component_moment_norm, component_orientations)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(out) :: &
            component_dipolar_moments
        logical, intent(in) :: dipolar
        class(Abstract_Component_Moment_Norm), intent(in) :: component_moment_norm
        class(Abstract_Component_Coordinates), intent(in) :: component_orientations

        if (dipolar) then
            allocate(Concrete_Component_Dipolar_Moments :: component_dipolar_moments)
        else
            allocate(Null_Component_Dipolar_Moments :: component_dipolar_moments)
        end if
        call component_dipolar_moments%construct(component_moment_norm, component_orientations)
    end subroutine allocate_and_construct_dipolar_moments

    subroutine allocate_and_construct_total_moment(component_total_moment, &
        component_dipolar_moments)
        class(Abstract_Component_Total_Moment), allocatable, intent(out) :: component_total_moment
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments

        if (component_is_dipolar(component_dipolar_moments)) then
            allocate(Concrete_Component_Total_Moment :: component_total_moment)
        else
            allocate(Null_Component_Total_Moment :: component_total_moment)
        end if
        call component_total_moment%construct(component_dipolar_moments)
    end subroutine allocate_and_construct_total_moment

    subroutine allocate_and_set_chemical_potential(component_chemical_potential, input_data, prefix)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(out) :: &
            component_chemical_potential
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, excess

        if (component_can_exchange(input_data, prefix)) then
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
    end subroutine allocate_and_set_chemical_potential

    subroutine component_factory_destroy_all(component)
        type(Component_Wrapper), intent(inout) :: component

        call component_factory_destroy(component%chemical_potential)
        call component_factory_destroy(component%total_moment)
        call component_factory_destroy(component%dipolar_moments)
        call component_factory_destroy(component%orientations)
        call component_factory_destroy(component%positions)
        call component_factory_destroy(component%moment_norm)
        call component_factory_destroy(component%wall_diameter)
        call component_factory_destroy(component%diameter)
        call component_factory_destroy(component%number)
    end subroutine component_factory_destroy_all

    subroutine deallocate_number(component_number)
        class(Abstract_Component_Number), allocatable, intent(inout) :: component_number

        if (allocated(component_number)) deallocate(component_number)
    end subroutine deallocate_number

    subroutine deallocate_diameter(component_diameter)
        class(Abstract_Component_Diameter), allocatable, intent(inout) :: component_diameter

        if (allocated(component_diameter)) deallocate(component_diameter)
    end subroutine deallocate_diameter

    subroutine deallocate_moment_norm(component_moment_norm)
        class(Abstract_Component_Moment_Norm), allocatable, intent(inout) :: component_moment_norm

        if (allocated(component_moment_norm)) deallocate(component_moment_norm)
    end subroutine deallocate_moment_norm

    subroutine destroy_and_deallocate_coordinates(component_coordinates)
        class(Abstract_Component_Coordinates), allocatable, intent(inout) :: component_coordinates

        call component_coordinates%destroy()
        if (allocated(component_coordinates)) deallocate(component_coordinates)
    end subroutine destroy_and_deallocate_coordinates

    subroutine destroy_and_deallocate_dipolar_moments(component_dipolar_moments)
        class(Abstract_Component_Dipolar_Moments), allocatable, intent(inout) :: &
            component_dipolar_moments

        call component_dipolar_moments%destroy()
        if (allocated(component_dipolar_moments)) deallocate(component_dipolar_moments)
    end subroutine destroy_and_deallocate_dipolar_moments

    subroutine destroy_and_deallocate_total_moment(component_total_moment)
        class(Abstract_Component_Total_Moment), allocatable, intent(inout) :: component_total_moment

        call component_total_moment%destroy()
        if (allocated(component_total_moment)) deallocate(component_total_moment)
    end subroutine destroy_and_deallocate_total_moment

    subroutine deallocate_chemical_potential(component_chemical_potential)
        class(Abstract_Component_Chemical_Potential), allocatable, intent(inout) :: &
            component_chemical_potential

        if (allocated(component_chemical_potential)) deallocate(component_chemical_potential)
    end subroutine deallocate_chemical_potential

end module procedures_component_factory
