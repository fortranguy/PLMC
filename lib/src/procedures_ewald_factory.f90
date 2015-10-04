module procedures_ewald_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use types_particles_wrapper, only: Particles_Wrapper
use procedures_property_inquirers, only: particles_are_dipolar
use types_potential_domain, only: Concrete_Potential_Domain
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, &
    Tabulated_Ewald_Real_Pair, Raw_Ewald_Real_Pair, Null_Ewald_Real_Pair
use class_ewald_real_particles, only: Abstract_Ewald_Real_Particles, &
    Concrete_Ewald_Real_Particles, Null_Ewald_Real_Particles
use types_ewald_wrapper, only: Ewald_Wrapper

implicit none

private
public :: ewald_factory_create, ewald_factory_destroy

interface ewald_factory_create
    module procedure :: ewald_factory_create_all
    module procedure :: allocate_and_construct_real_pair
    module procedure :: allocate_and_construct_real_particles
end interface ewald_factory_create

interface ewald_factory_destroy
    module procedure :: destroy_and_deallocate_real_particles
    module procedure :: destroy_and_deallocate_real_pair
    module procedure :: ewald_factory_destroy_all
end interface ewald_factory_destroy

contains

    subroutine ewald_factory_create_all(ewald, input_data, prefix, environment, particles)
        type(Ewald_Wrapper), intent(out) :: ewald
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        type(Environment_Wrapper), intent(in) :: environment
        type(Particles_Wrapper), intent(in) :: particles

        real(DP) :: alpha

        call set_alpha(alpha, input_data, prefix, environment%periodic_box)
        call ewald_factory_create(ewald%real_pair, input_data, prefix//"Real.", alpha, &
            environment%periodic_box, particles%diameter, particles%dipolar_moments)
        call ewald_factory_create(ewald%real_particles, environment%periodic_box, &
            particles%positions, particles%dipolar_moments)
    end subroutine ewald_factory_create_all

    subroutine set_alpha(alpha, input_data, prefix, periodic_box)
        real(DP), intent(out) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Periodic_Box), intent(in) :: periodic_box

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: box_size(num_dimensions), alpha_times_box

        data_field = prefix//"alpha * box size(1)"
        call input_data%get(data_field, alpha_times_box, data_found)
        call check_data_found(data_field, data_found)
        box_size = periodic_box%get_size()
        alpha = alpha_times_box / box_size(1)
        deallocate(data_field)
    end subroutine set_alpha

    subroutine allocate_and_construct_real_pair(real_pair, input_data, prefix, alpha, &
        periodic_box, particles_diameter, particles_dipolar_moments)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        call allocate_real_pair(real_pair, input_data, prefix, particles_dipolar_moments)
        call construct_real_pair(real_pair, input_data, prefix, alpha, periodic_box, &
            particles_diameter, particles_dipolar_moments)
    end subroutine allocate_and_construct_real_pair

    subroutine allocate_real_pair(real_pair, input_data, prefix, particles_dipolar_moments)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (particles_are_dipolar(particles_dipolar_moments)) then
            data_field = prefix//"tabulated"
            call input_data%get(data_field, tabulated_potential, data_found)
            call check_data_found(data_field, data_found)
            if(tabulated_potential) then
                allocate(Tabulated_Ewald_Real_Pair :: real_pair)
            else
                allocate(Raw_Ewald_Real_Pair :: real_pair)
            end if
            deallocate(data_field)
        else
            allocate(Null_Ewald_Real_Pair :: real_pair)
        end if
    end subroutine allocate_real_pair

    subroutine construct_real_pair(real_pair, input_data, prefix, alpha, periodic_box, &
        particles_diameter, particles_dipolar_moments)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: real_pair
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: box_size(num_dimensions), max_over_box

        if (particles_are_dipolar(particles_dipolar_moments)) then
            domain%min = particles_diameter%get_min()
            data_field = prefix//"max distance / box size(1)"
            call input_data%get(data_field, max_over_box, data_found)
            call check_data_found(data_field, data_found)
            box_size = periodic_box%get_size()
            domain%max = max_over_box * box_size(1)
            select type (real_pair)
                type is (Tabulated_Ewald_Real_Pair)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call real_pair%construct(domain, alpha)
    end subroutine construct_real_pair

    subroutine allocate_and_construct_real_particles(real_particles, periodic_box, &
        particles_positions, particles_dipolar_moments)
        class(Abstract_Ewald_Real_Particles), allocatable, intent(out) :: real_particles
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments

        if (particles_are_dipolar(particles_dipolar_moments)) then
            allocate(Concrete_Ewald_Real_Particles :: real_particles)
        else
            allocate(Null_Ewald_Real_Particles :: real_particles)
        end if
        call real_particles%construct(periodic_box, particles_positions, particles_dipolar_moments)
    end subroutine allocate_and_construct_real_particles

    subroutine ewald_factory_destroy_all(ewald)
        type(Ewald_Wrapper), intent(inout) :: ewald

        call ewald_factory_destroy(ewald%real_particles)
        call ewald_factory_destroy(ewald%real_pair)
    end subroutine ewald_factory_destroy_all

    subroutine destroy_and_deallocate_real_particles(real_particles)
        class(Abstract_Ewald_Real_Particles), allocatable, intent(inout) :: real_particles

        call real_particles%destroy()
        if (allocated(real_particles)) deallocate(real_particles)
    end subroutine destroy_and_deallocate_real_particles

    subroutine destroy_and_deallocate_real_pair(real_pair)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(inout) :: real_pair

        call real_pair%destroy()
        if (allocated(real_pair)) deallocate(real_pair)
    end subroutine destroy_and_deallocate_real_pair

end module procedures_ewald_factory
