module procedures_short_potential_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_errors, only: error_exit
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_diameter, only: Abstract_Particles_Diameter, &
    Null_Particles_Diameter
use procedures_particles_factory, only: particles_exist
use class_potential_expression, only: Abstract_Potential_Expression, &
    Null_Potential_Expression, Lennard_Jones_Expression
use types_potential_domain, only: Concrete_Potential_Domain
use class_pair_potential, only: Abstract_Pair_Potential, &
    Null_Pair_Potential, Tabulated_Pair_Potential, Raw_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential, &
    Concrete_Particles_Potential, Null_Particles_Potential
use types_short_potential, only: Short_Potential_Wrapper

implicit none

private
public :: allocate_and_set_expression, allocate_and_construct_pair, &
    allocate_and_construct_particles_potential

contains

    subroutine potential_factory_construct(short_potential, input_data, prefix, periodic_box, &
        particles_diameter)
        type(Short_Potential_Wrapper), intent(out) :: short_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        call allocate_and_set_expression(short_potential%expression, input_data, prefix, &
            particles_diameter)
        call allocate_and_construct_pair(short_potential%pair, input_data, prefix, &
            particles_diameter, short_potential%expression)
        call allocate_and_construct_particles_potential(short_potential%particles, &
            periodic_box, particles_diameter)
    end subroutine potential_factory_construct

    subroutine allocate_and_set_expression(potential_expression, input_data, prefix, &
        particles_diameter)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        call allocate_expression(potential_expression, input_data, prefix, particles_diameter)
        call set_expression(potential_expression, input_data, prefix)
    end subroutine allocate_and_set_expression

    subroutine allocate_expression(potential_expression, input_data, prefix, particles_diameter)
        class(Abstract_Potential_Expression), allocatable, intent(out) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        character(len=:), allocatable :: data_field, potential_name
        logical :: data_found

        if (particles_exist(particles_diameter)) then
            data_field = prefix//".Potential.name"
            call input_data%get(data_field, potential_name, data_found)
            call test_data_found(data_field, data_found)
            select case(potential_name)
                case ("null")
                    allocate(Null_Potential_Expression :: potential_expression)
                case ("LJ")
                    allocate(Lennard_Jones_Expression :: potential_expression)
                case default
                    call error_exit(potential_name//" unknown potential_name."//&
                        "Choose between: 'null' and LJ.")
            end select
            deallocate(potential_name)
            deallocate(data_field)
        else
            allocate(Null_Potential_Expression :: potential_expression)
        end if
    end subroutine allocate_expression

    subroutine set_expression(potential_expression, input_data, prefix)
        class(Abstract_Potential_Expression), intent(inout) :: potential_expression
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: LJ_epsilon, LJ_sigma

        select type(potential_expression)
            type is (Null_Potential_Expression)
                call potential_expression%set()
            type is (Lennard_Jones_Expression)
                data_field = prefix//".Potential.epsilon"
                call input_data%get(data_field, LJ_epsilon, data_found)
                call test_data_found(data_field, data_found)
                data_field = prefix//".Potential.sigma"
                call input_data%get(data_field, LJ_sigma, data_found)
                call test_data_found(data_field, data_found)
                call potential_expression%set(LJ_epsilon, LJ_sigma)
        end select
        if (allocated(data_field)) deallocate(data_field)
    end subroutine set_expression

    subroutine allocate_and_construct_pair(pair_potential, input_data, prefix, &
        particles_diameter, potential_expression)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression

        call allocate_pair(pair_potential, input_data, prefix, particles_diameter)
        call construct_pair(pair_potential, input_data, prefix, particles_diameter, &
            potential_expression)
    end subroutine allocate_and_construct_pair

    subroutine allocate_pair(pair_potential, input_data, prefix, particles_diameter)
        class(Abstract_Pair_Potential), allocatable, intent(out) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (particles_exist(particles_diameter)) then
            data_field = prefix//".Potential.tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call test_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_Pair_Potential :: pair_potential)
            else
                allocate(Raw_Pair_Potential :: pair_potential)
            end if
            deallocate(data_field)
        else
            allocate(Null_Pair_Potential :: pair_potential)
        end if
    end subroutine allocate_pair

    subroutine construct_pair(pair_potential, input_data, prefix, particles_diameter, &
        potential_expression)
        class(Abstract_Pair_Potential), intent(inout) :: pair_potential
        type(json_file), target, intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Potential_Expression), intent(in) :: potential_expression

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain

        if (particles_exist(particles_diameter)) then
            domain%min = particles_diameter%get_min()
            data_field = prefix//".Potential.max distance"
            call input_data%get(data_field, domain%max, data_found)
            call test_data_found(data_field, data_found)
            select type (pair_potential)
                type is (Tabulated_Pair_Potential)
                    data_field = prefix//".Potential.delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call test_data_found(data_field, data_found)
            end select
        end if
        call pair_potential%construct(domain, potential_expression)
    end subroutine construct_pair

    subroutine allocate_and_construct_particles_potential(particles_potential, &
        periodic_box, particles_diameter)
        class(Abstract_Particles_Potential), allocatable, intent(out) :: particles_potential
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter

        if (particles_exist(particles_diameter)) then
            allocate(Concrete_Particles_Potential :: particles_potential)
        else
            allocate(Null_Particles_Potential :: particles_potential)
        end if
        call particles_potential%construct(periodic_box)
    end subroutine allocate_and_construct_particles_potential

    subroutine potential_factory_destroy(short_potential)
        type(Short_Potential_Wrapper), intent(inout) :: short_potential

        call short_potential%particles%destroy()
        if (allocated(short_potential%particles)) deallocate(short_potential%particles)
        call short_potential%pair%destroy()
        if (allocated(short_potential%pair)) deallocate(short_potential%pair)
        if (allocated(short_potential%expression)) deallocate(short_potential%expression)
    end subroutine potential_factory_destroy

end module procedures_short_potential_factory
