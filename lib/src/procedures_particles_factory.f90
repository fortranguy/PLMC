module procedures_particles_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_errors, only: error_exit
use procedures_coordinates, only: read_coordinates
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_number, only: Abstract_Particles_Number, &
    Concrete_Particles_Number, Null_Particles_Number
use class_particles_diameter, only: Abstract_Particles_Diameter, &
    Concrete_Particles_Diameter, Null_Particles_Diameter
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm, &
    Concrete_Particles_Moment_Norm, Null_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions, &
    Concrete_Particles_Positions, Null_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations, &
    Concrete_Particles_Orientations, Null_Particles_Orientations
use class_particles_chemical_potential, only : Abstract_Particles_Chemical_Potential, &
    Concrete_Particles_Chemical_Potential, Null_Particles_Chemical_Potential
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments, &
    Concrete_Particles_Dipolar_Moments, Null_Particles_Dipolar_Moments
use class_particles_total_moment, only: Abstract_Particles_Total_Moment, &
    Concrete_Particles_Total_Moment, Null_Particles_Total_Moment
use types_particles, only: Particles_Wrapper, Mixture_Wrapper
use procedures_types_selectors, only: particles_exist, particles_have_positions, &
    particles_have_orientations, particles_are_dipolar

implicit none

private
public :: mixture_factory_create, mixture_factory_destroy, &
    particles_factory_create, particles_factory_destroy, &
    allocate_and_set_number, allocate_and_set_diameter, allocate_and_set_inter_diameter, &
    allocate_and_construct_positions, set_positions

contains

    subroutine mixture_factory_create(mixture, input_data, prefix, periodic_box)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call particles_factory_create(mixture%components(1), input_data, prefix//".Component 1", &
            periodic_box)
        call particles_factory_create(mixture%components(2), input_data, prefix//".Component 2", &
            periodic_box)
        call allocate_and_set_inter_diameter(mixture%inter_diameters, &
            mixture%components(1)%diameter, mixture%components(2)%diameter, input_data, &
            prefix//".Inter 12")
    end subroutine mixture_factory_create

    subroutine particles_factory_create(particles, input_data, prefix, periodic_box)
        type(Particles_Wrapper), intent(out) :: particles
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        call allocate_and_set_number(particles%number, input_data, prefix)
        call allocate_and_set_diameter(particles%diameter, input_data, prefix)
        call allocate_and_set_moment_norm(particles%moment_norm, input_data, prefix)
        call allocate_and_construct_positions(particles%positions, periodic_box, particles%number)
        call set_positions(particles%positions, input_data, prefix)
        call allocate_and_construct_orientations(particles%orientations, input_data, prefix, &
            particles%number)
        call set_orientations(particles%orientations, input_data, prefix)
        call allocate_and_construct_dipolar_moments(particles%dipolar_moments, &
            particles%moment_norm, particles%orientations)
        call allocate_and_construct_total_moment(particles%total_moment, particles%dipolar_moments)
        call allocate_and_set_chemical_potential(particles%chemical_potential, input_data, prefix)
    end subroutine particles_factory_create

    subroutine allocate_and_set_number(particles_number, input_data, prefix)
        class(Abstract_Particles_Number), allocatable, intent(out) :: particles_number
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: num_particles

        if (particles_exist_from_json(input_data, prefix)) then
            data_field = prefix//".number"
            call input_data%get(data_field, num_particles, data_found)
            call test_data_found(data_field, data_found)
            allocate(Concrete_Particles_Number :: particles_number)
            deallocate(data_field)
        else
            allocate(Null_Particles_Number :: particles_number)
        end if
        call particles_number%set(num_particles)
    end subroutine allocate_and_set_number

    logical function particles_exist_from_json(input_data, prefix) result(particles_exist)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        data_field = prefix//".exist"
        call input_data%get(data_field, particles_exist, data_found)
        call test_data_found(data_field, data_found)
        deallocate(data_field)
    end function particles_exist_from_json

    subroutine allocate_and_set_diameter(particles_diameter, input_data, prefix)
        class(Abstract_Particles_Diameter), allocatable, intent(out) :: particles_diameter
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: diameter, diameter_min_factor

        if (particles_exist_from_json(input_data, prefix)) then
            data_field = prefix//".diameter"
            call input_data%get(data_field, diameter, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".minimum diameter factor"
            call input_data%get(data_field, diameter_min_factor, data_found)
            call test_data_found(data_field, data_found)
            allocate(Concrete_Particles_Diameter :: particles_diameter)
            deallocate(data_field)
        else
            allocate(Null_Particles_Diameter :: particles_diameter)
        end if
        call particles_diameter%set(diameter, diameter_min_factor)
    end subroutine allocate_and_set_diameter

    subroutine allocate_and_set_inter_diameter(inter_particles_diameter, &
        particles_diameter_1, particles_diameter_2, input_data, prefix)
        class(Abstract_Particles_Diameter), allocatable, intent(out) :: inter_particles_diameter
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter_1, particles_diameter_2
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: inter_diameter, inter_diameter_min_factor
        real(DP) :: inter_diameter_offset

        if (particles_exist(particles_diameter_1) .and. particles_exist(particles_diameter_2)) then
            data_field = prefix//".offset"
            call input_data%get(data_field, inter_diameter_offset, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".minimum diameter factor"
            call input_data%get(data_field, inter_diameter_min_factor, data_found)
            call test_data_found(data_field, data_found)
            allocate(Concrete_Particles_Diameter :: inter_particles_diameter)
            deallocate(data_field)
        else
            allocate(Null_Particles_Diameter :: inter_particles_diameter)
        end if
        inter_diameter = (particles_diameter_1%get() + particles_diameter_2%get()) / 2._DP + &
            inter_diameter_offset
        call inter_particles_diameter%set(inter_diameter, inter_diameter_min_factor)
    end subroutine allocate_and_set_inter_diameter

    subroutine allocate_and_set_moment_norm(particles_moment_norm, input_data, prefix)
        class(Abstract_Particles_Moment_Norm), allocatable, intent(out) :: particles_moment_norm
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: moment_norm

        if (particles_are_dipolar_from_json(input_data, prefix)) then
            data_field = prefix//".moment norm"
            call input_data%get(data_field, moment_norm, data_found)
            call test_data_found(data_field, data_found)
            allocate(Concrete_Particles_Moment_Norm :: particles_moment_norm)
            deallocate(data_field)
        else
            allocate(Null_Particles_Moment_Norm :: particles_moment_norm)
        end if
        call particles_moment_norm%set(moment_norm)
    end subroutine allocate_and_set_moment_norm

    logical function particles_are_dipolar_from_json(input_data, prefix) &
        result(particles_are_dipolar)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (particles_exist_from_json(input_data, prefix)) then
            data_field = prefix//".are dipolar"
            call input_data%get(data_field, particles_are_dipolar, data_found)
            call test_data_found(data_field, data_found)
            deallocate(data_field)
        else
            particles_are_dipolar = .false.
        end if
    end function particles_are_dipolar_from_json

    subroutine allocate_and_construct_positions(particles_positions, periodic_box, particles_number)
        class(Abstract_Particles_Positions), allocatable, intent(out) :: particles_positions
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Number), intent(in) :: particles_number

        if (particles_exist(particles_number)) then
            allocate(Concrete_Particles_Positions :: particles_positions)
        else
            allocate(Null_Particles_Positions :: particles_positions)
        end if
        call particles_positions%construct(periodic_box, particles_number)
    end subroutine allocate_and_construct_positions

    subroutine allocate_and_construct_orientations(particles_orientations, input_data, prefix, &
        particles_number)
        class(Abstract_Particles_Orientations), allocatable, intent(out) :: particles_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Number), intent(in) :: particles_number

        if (particles_are_dipolar_from_json(input_data, prefix)) then
            allocate(Concrete_Particles_Orientations :: particles_orientations)
        else
            allocate(Null_Particles_Orientations :: particles_orientations)
        end if
        call particles_orientations%construct(particles_number)
    end subroutine allocate_and_construct_orientations

    subroutine set_positions(particles_positions, input_data, prefix)
        class(Abstract_Particles_Positions), intent(inout) :: particles_positions
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, filename
        logical :: data_found
        real(DP), allocatable :: file_positions(:, :)
        integer :: i_particle

        if (.not. particles_have_positions(particles_positions)) return
        data_field = prefix//".initial positions"
        call input_data%get(data_field, filename, data_found)
        call test_data_found(data_field, data_found)
        call read_coordinates(file_positions, filename)
        if (size(file_positions, 2) /= particles_positions%get_num()) then
            call error_exit("set_positions from "//filename//": wrong number of lines.")
        end if
        do i_particle = 1, particles_positions%get_num()
            call particles_positions%set(i_particle, file_positions(:, i_particle))
        end do
        if (allocated(file_positions)) deallocate(file_positions)
        deallocate(filename)
        deallocate(data_field)
    end subroutine set_positions

    subroutine set_orientations(particles_orientations, input_data, prefix)
        class(Abstract_Particles_Orientations), intent(inout) :: particles_orientations
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field, filename
        logical :: data_found
        real(DP), allocatable :: file_orientations(:, :)
        integer :: i_particle

        if (.not. particles_have_orientations(particles_orientations)) return
        data_field = prefix//".initial orientations"
        call input_data%get(data_field, filename, data_found)
        call test_data_found(data_field, data_found)
        call read_coordinates(file_orientations, filename)
        if (size(file_orientations, 2) /= particles_orientations%get_num()) then
            call error_exit("set_orientations from "//filename//": wrong number of lines.")
        end if
        do i_particle = 1, particles_orientations%get_num()
            call particles_orientations%set(i_particle, file_orientations(:, i_particle))
        end do
        if (allocated(file_orientations)) deallocate(file_orientations)
        deallocate(filename)
        deallocate(data_field)
    end subroutine set_orientations

    subroutine allocate_and_construct_dipolar_moments(particles_dipolar_moments, &
        particles_moment_norm, particles_orientations)
        class(Abstract_Particles_Dipolar_Moments), allocatable, intent(out) :: &
            particles_dipolar_moments
        class(Abstract_Particles_Moment_Norm), intent(in) :: particles_moment_norm
        class(Abstract_Particles_Orientations), intent(in) :: particles_orientations

        if (particles_are_dipolar(particles_moment_norm, particles_orientations)) then
            allocate(Concrete_Particles_Dipolar_Moments :: particles_dipolar_moments)
        else
            allocate(Null_Particles_Dipolar_Moments :: particles_dipolar_moments)
        end if
        call particles_dipolar_moments%construct(particles_moment_norm, particles_orientations)
    end subroutine allocate_and_construct_dipolar_moments

    subroutine allocate_and_construct_total_moment(particles_total_moment, &
        particles_dipolar_moments)
        class(Abstract_Particles_Total_Moment), allocatable, intent(out) :: particles_total_moment
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        if (particles_are_dipolar(particles_dipolar_moments)) then
            allocate(Concrete_Particles_Total_Moment :: particles_total_moment)
        else
            allocate(Null_Particles_Total_Moment :: particles_total_moment)
        end if
        call particles_total_moment%construct(particles_dipolar_moments)
    end subroutine allocate_and_construct_total_moment

    subroutine allocate_and_set_chemical_potential(particles_chemical_potential, input_data, prefix)
        class(Abstract_Particles_Chemical_Potential), allocatable, intent(out) :: &
            particles_chemical_potential
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: density, excess

        if (particles_can_exchange_from_json(input_data, prefix)) then
            data_field = prefix//".Chemical Potential.density"
            call input_data%get(data_field, density, data_found)
            call test_data_found(data_field, data_found)
            data_field = prefix//".Chemical Potential.excess"
            call input_data%get(data_field, excess, data_found)
            call test_data_found(data_field, data_found)
            allocate(Concrete_Particles_Chemical_Potential :: particles_chemical_potential)
            deallocate(data_field)
        else
            allocate(Null_Particles_Chemical_Potential :: particles_chemical_potential)
        end if
        call particles_chemical_potential%set(density, excess)
    end subroutine allocate_and_set_chemical_potential

    logical function particles_can_exchange_from_json(input_data, prefix) &
        result(particles_exchange_with_reservoir)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found

        if (particles_exist_from_json(input_data, prefix)) then
            data_field = prefix//".can exchange"
            call input_data%get(data_field, particles_exchange_with_reservoir, data_found)
            call test_data_found(data_field, data_found)
            deallocate(data_field)
        else
            particles_exchange_with_reservoir = .false.
        end if
    end function particles_can_exchange_from_json

    subroutine mixture_factory_destroy(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        if (allocated(mixture%inter_diameters)) deallocate(mixture%inter_diameters)
        call particles_factory_destroy(mixture%components(2))
        call particles_factory_destroy(mixture%components(1))
    end subroutine mixture_factory_destroy

    subroutine particles_factory_destroy(particles)
        type(Particles_Wrapper), intent(inout) :: particles

        if (allocated(particles%chemical_potential)) deallocate(particles%chemical_potential)
        call particles%total_moment%destroy()
        if (allocated(particles%total_moment)) deallocate(particles%total_moment)
        call particles%dipolar_moments%destroy()
        if (allocated(particles%dipolar_moments)) deallocate(particles%dipolar_moments)
        call particles%orientations%destroy()
        if (allocated(particles%orientations)) deallocate(particles%orientations)
        call particles%positions%destroy()
        if (allocated(particles%positions)) deallocate(particles%positions)
        if (allocated(particles%moment_norm)) deallocate(particles%moment_norm)
        if (allocated(particles%diameter)) deallocate(particles%diameter)
        if (allocated(particles%number)) deallocate(particles%number)
    end subroutine particles_factory_destroy

end module procedures_particles_factory
