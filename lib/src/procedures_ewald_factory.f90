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
use types_particles_wrapper, only: Particles_Wrapper, Mixture_Wrapper
use types_potential_domain, only: Concrete_Potential_Domain
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, &
    Tabulated_Ewald_Real_Pair, Raw_Ewald_Real_Pair, Null_Ewald_Real_Pair
use class_ewald_real_particles, only: Abstract_Ewald_Real_Particles, &
    Concrete_Ewald_Real_Particles, Null_Ewald_Real_Particles
use types_ewald_wrapper, only: Ewald_Wrapper, Ewald_Wrapper_Macro, Ewald_Wrapper_Micro
use procedures_property_inquirers, only: particles_are_dipolar, particles_interact

implicit none

private
public :: ewald_factory_create, ewald_factory_destroy

interface ewald_factory_create
    module procedure :: ewald_factory_create_all
    module procedure :: ewald_factory_create_macro
    module procedure :: ewald_factory_create_micro
    module procedure :: allocate_and_construct_real_pair
    module procedure :: allocate_and_construct_real_particles
end interface ewald_factory_create

interface ewald_factory_destroy
    module procedure :: destroy_and_deallocate_real_particles
    module procedure :: destroy_and_deallocate_real_pair
    module procedure :: ewald_factory_destroy_micro
    module procedure :: ewald_factory_destroy_macro
    module procedure :: ewald_factory_destroy_all
end interface ewald_factory_destroy

contains

    subroutine ewald_factory_create_all(ewald, environment, particles, input_data, prefix)
        type(Ewald_Wrapper), intent(out) :: ewald
        type(Environment_Wrapper), intent(in) :: environment
        type(Particles_Wrapper), intent(in) :: particles
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: alpha
        logical :: dipolar

        dipolar = particles_are_dipolar(particles%dipolar_moments)
        call set_alpha(alpha, dipolar, environment%periodic_box, input_data, prefix)
        call ewald_factory_create(ewald%real_pair, dipolar, alpha, environment%periodic_box, &
            particles%diameter, input_data, prefix//"Real.")
        call ewald_factory_create(ewald%real_particles, dipolar, environment%periodic_box, &
            particles%positions, particles%dipolar_moments, ewald%real_pair)
    end subroutine ewald_factory_create_all

    subroutine ewald_factory_create_macro(ewald_macro, ewald_micro, environment, particles)
        type(Ewald_Wrapper_Macro), intent(out) :: ewald_macro
        type(Ewald_Wrapper_Micro), intent(in) :: ewald_micro
        type(Environment_Wrapper), intent(in) :: environment
        type(Particles_Wrapper), intent(in) :: particles

        logical :: interact

        interact = particles_interact(ewald_micro%real_pair)
        call ewald_factory_create(ewald_macro%real_particles, interact, environment%periodic_box, &
            particles%positions, particles%dipolar_moments, ewald_micro%real_pair)
    end subroutine ewald_factory_create_macro

    subroutine ewald_factory_create_micro(ewald_micro, dipolar, environment, particles_diameter, &
        input_data, prefix)
        type(Ewald_Wrapper_Micro), intent(out) :: ewald_micro
        logical, intent(in) :: dipolar
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: alpha

        call set_alpha(alpha, dipolar, environment%periodic_box, input_data, prefix)
        call ewald_factory_create(ewald_micro%real_pair, dipolar, alpha, &
            environment%periodic_box, particles_diameter, input_data, prefix//"Real.")
    end subroutine ewald_factory_create_micro

    subroutine set_alpha(alpha, dipolar, periodic_box, input_data, prefix)
        real(DP), intent(out) :: alpha
        logical, intent(in) :: dipolar
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: box_size(num_dimensions), alpha_times_box

        if (dipolar) then
            data_field = prefix//"alpha times box edge"
            call input_data%get(data_field, alpha_times_box, data_found)
            call check_data_found(data_field, data_found)
            box_size = periodic_box%get_size()
            alpha = alpha_times_box / box_size(1)
            deallocate(data_field)
        else
            alpha = 0._DP
        end if
    end subroutine set_alpha

    subroutine allocate_and_construct_real_pair(real_pair, dipolar, alpha, periodic_box, &
        particles_diameter, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: dipolar
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_real_pair(real_pair, dipolar, input_data, prefix)
        call construct_real_pair(real_pair, dipolar, alpha, periodic_box, particles_diameter, &
            input_data, prefix)
    end subroutine allocate_and_construct_real_pair

    subroutine allocate_real_pair(real_pair, dipolar, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: dipolar
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (dipolar) then
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

    subroutine construct_real_pair(real_pair, dipolar, alpha, periodic_box, particles_diameter, &
        input_data, prefix)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: real_pair
        logical, intent(in) :: dipolar
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Diameter), intent(in) :: particles_diameter
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: box_size(num_dimensions), max_over_box

        if (dipolar) then
            domain%min = particles_diameter%get_min()
            data_field = prefix//"max distance over box edge"
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

    subroutine allocate_and_construct_real_particles(real_particles, dipolar, periodic_box, &
        particles_positions, particles_dipolar_moments, real_pair)
        class(Abstract_Ewald_Real_Particles), allocatable, intent(out) :: real_particles
        logical, intent(in) :: dipolar
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Particles_Positions), intent(in) :: particles_positions
        class(Abstract_Particles_Dipolar_Moments), intent(in) :: particles_dipolar_moments
        class(Abstract_Ewald_Real_Pair), intent(in) :: real_pair

        if (dipolar) then
            allocate(Concrete_Ewald_Real_Particles :: real_particles)
        else
            allocate(Null_Ewald_Real_Particles :: real_particles)
        end if
        call real_particles%construct(periodic_box, particles_positions, &
            particles_dipolar_moments, real_pair)
    end subroutine allocate_and_construct_real_particles

    subroutine ewald_factory_destroy_all(ewald)
        type(Ewald_Wrapper), intent(inout) :: ewald

        call ewald_factory_destroy(ewald%real_particles)
        call ewald_factory_destroy(ewald%real_pair)
    end subroutine ewald_factory_destroy_all

    subroutine ewald_factory_destroy_micro(ewald_micro)
        type(Ewald_Wrapper_Micro), intent(inout) :: ewald_micro

        call ewald_factory_destroy(ewald_micro%real_pair)
    end subroutine ewald_factory_destroy_micro

    subroutine ewald_factory_destroy_macro(ewald_macro)
        type(Ewald_Wrapper_Macro), intent(inout) :: ewald_macro

        call ewald_factory_destroy(ewald_macro%real_particles)
    end subroutine ewald_factory_destroy_macro

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
