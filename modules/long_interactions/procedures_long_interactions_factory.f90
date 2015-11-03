module procedures_long_interactions_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_errors, only: warning_continue
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_environment_wrapper, only: Environment_Wrapper
use class_minimum_distance, only: Abstract_Minimum_Distance
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Minimum_Distances_Wrapper, Mixture_Wrapper
use types_potential_domain, only: Concrete_Potential_Domain
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, Tabulated_Ewald_Real_Pair, &
    Raw_Ewald_Real_Pair, Null_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component, &
    Concrete_Ewald_Real_Component, Null_Ewald_Real_Component
use class_ewald_reci_structures, only: Abstract_Ewald_Reci_Structures
use class_ewald_real_visitor, only: Abstract_Ewald_Real_Visitor, Concrete_Ewald_Real_Visitor, &
    Null_Ewald_Real_Visitor
use types_long_interactions_wrapper, only: Ewald_Real_Pair_Wrapper, Ewald_Real_Pairs_Wrapper, &
    Ewald_Real_Component_Wrapper, Long_Interactions_Wrapper
use procedures_property_inquirers, only: use_permittivity, use_reciprocal_lattice, &
    component_is_dipolar, components_interact

implicit none

private
public :: long_interactions_create, long_interactions_set, long_interactions_destroy

interface long_interactions_create
    module procedure :: create_all
    module procedure :: create_real_visitor
    module procedure :: create_real_components
    module procedure :: create_real_component
    module procedure :: create_real_pairs
    module procedure :: create_real_pair
end interface long_interactions_create

interface long_interactions_set
    module procedure :: set_alpha
    module procedure :: set_are_dipolar
end interface long_interactions_set

interface long_interactions_check
    module procedure :: check_consistency
end interface

interface long_interactions_destroy
    module procedure :: destroy_real_pair
    module procedure :: destroy_real_pairs
    module procedure :: destroy_real_component
    module procedure :: destroy_real_components
    module procedure :: destroy_real_visitor
    module procedure :: destroy_all
end interface long_interactions_destroy

contains

    subroutine create_all(long_interactions, environment, mixture, input_data, prefix)
        type(Long_Interactions_Wrapper), intent(out) :: long_interactions
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: are_dipolar(size(mixture%components))
        real(DP) :: alpha

        call long_interactions_set(are_dipolar, mixture%components)
        call long_interactions_check(any(are_dipolar), environment%reciprocal_lattice, &
            environment%permittivity)
        call long_interactions_create(long_interactions%real_visitor, environment%periodic_box, &
            any(are_dipolar))
        call long_interactions_set(alpha, environment%periodic_box, any(are_dipolar), input_data, &
            prefix)
        call long_interactions_create(long_interactions%real_pairs, environment, mixture%&
            components_min_distances, are_dipolar, alpha, input_data, prefix//"Real.")
        call long_interactions_create(long_interactions%real_components, environment%periodic_box, &
            mixture%components, long_interactions%real_pairs)
    end subroutine create_all

    subroutine destroy_all(long_interactions)
        type(Long_Interactions_Wrapper), intent(inout) :: long_interactions

        call long_interactions_destroy(long_interactions%real_pairs)
        call long_interactions_destroy(long_interactions%real_components)
        call long_interactions_destroy(long_interactions%real_visitor)
    end subroutine destroy_all

    subroutine create_real_visitor(visitor, periodic_box, dipoles_exist)
        class(Abstract_Ewald_Real_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            allocate(Concrete_Ewald_Real_Visitor :: visitor)
        else
            allocate(Null_Ewald_Real_Visitor :: visitor)
        end if
        call visitor%construct(periodic_box)
    end subroutine create_real_visitor

    subroutine destroy_real_visitor(visitor)
        class(Abstract_Ewald_Real_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_real_visitor

    subroutine create_real_components(real_components, periodic_box, components, real_pairs)
        type(Ewald_Real_Component_Wrapper), allocatable, intent(out) :: real_components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)
        type(Ewald_Real_Pairs_Wrapper), intent(in) :: real_pairs(:)

        integer :: j_component, i_component
        integer :: j_pair, i_pair

        allocate(real_components(size(real_pairs), size(real_pairs)))

        do j_component = 1, size(real_components, 2)
            do i_component = 1, size(real_components, 1)
                j_pair = maxval([j_component, i_component])
                i_pair = minval([j_component, i_component])
                associate (pair_ij => real_pairs(j_pair)%with_components(i_pair)%real_pair)
                    call long_interactions_create(real_components(i_component, j_component)%&
                        real_component, periodic_box, components(i_component)%positions, &
                        components(i_component)%dipolar_moments, pair_ij)
                end associate
            end do
        end do
    end subroutine create_real_components

    subroutine destroy_real_components(real_components)
        type(Ewald_Real_Component_Wrapper), allocatable, intent(inout) :: real_components(:, :)

        integer :: j_component, i_component

        if (allocated(real_components)) then
            do j_component = size(real_components, 2), 1, -1
                do i_component = size(real_components, 1), 1, -1
                    call long_interactions_destroy(real_components(i_component, j_component)%&
                        real_component)
                end do
            end do
            deallocate(real_components)
        end if
    end subroutine destroy_real_components

    subroutine create_real_component(real_component, periodic_box, positions, dipolar_moments, &
        real_pair)
        class(Abstract_Ewald_Real_Component), allocatable, intent(out) :: real_component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        class(Abstract_Ewald_Real_Pair), intent(in) :: real_pair

        if (components_interact(real_pair)) then
            allocate(Concrete_Ewald_Real_Component :: real_component)
        else
            allocate(Null_Ewald_Real_Component :: real_component)
        end if
        call real_component%construct(periodic_box, positions, dipolar_moments, real_pair)
    end subroutine create_real_component

    subroutine destroy_real_component(real_component)
        class(Abstract_Ewald_Real_Component), allocatable, intent(inout) :: real_component

        if (allocated(real_component)) then
            call real_component%destroy()
            deallocate(real_component)
        end if
    end subroutine destroy_real_component

    subroutine create_real_pairs(real_pairs, environment, min_distances, are_dipolar, &
        alpha, input_data, prefix)
        type(Ewald_Real_Pairs_Wrapper), allocatable, intent(out) :: real_pairs(:)
        type(Environment_Wrapper), intent(in) :: environment
        type(Minimum_Distances_Wrapper), intent(in) :: min_distances(:)
        logical, intent(in) :: are_dipolar(:)
        real(DP), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: interact_ij
        integer :: j_component, i_component

        allocate(real_pairs(size(are_dipolar)))
        do j_component = 1, size(are_dipolar)
            allocate(real_pairs(j_component)%with_components(j_component))
            do i_component = 1, size(real_pairs(j_component)%with_components)
                interact_ij = are_dipolar(i_component) .and. are_dipolar(j_component)
                associate (min_distance_ij => min_distances(j_component)%&
                    with_components(i_component)%min_distance)
                    call long_interactions_create(real_pairs(j_component)%&
                        with_components(i_component)%real_pair, environment, min_distance_ij, &
                        interact_ij, alpha, input_data, prefix)
                end associate
            end do
        end do
    end subroutine create_real_pairs

    subroutine destroy_real_pairs(real_pairs)
        type(Ewald_Real_Pairs_Wrapper), allocatable, intent(inout) :: real_pairs(:)

        integer :: j_component, i_component

        if (allocated(real_pairs)) then
            do j_component = size(real_pairs), 1, -1
                if (allocated(real_pairs(j_component)%with_components)) then
                    do i_component = size(real_pairs(j_component)%with_components), 1, -1
                        call long_interactions_destroy(real_pairs(j_component)%&
                            with_components(i_component)%real_pair)
                    end do
                    deallocate(real_pairs(j_component)%with_components)
                end if
            end do
            deallocate(real_pairs)
        end if
    end subroutine destroy_real_pairs

    subroutine create_real_pair(real_pair, environment, min_distance, interact, alpha, &
        input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        real(DP), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_real_pair(real_pair, interact, input_data, prefix)
        call construct_real_pair(real_pair, environment, min_distance, interact, alpha, &
            input_data, prefix)
    end subroutine create_real_pair

    subroutine allocate_real_pair(real_pair, interact, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: interact
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (interact) then
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

    subroutine construct_real_pair(real_pair, environment, min_distance, interact, alpha, &
        input_data, prefix)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: real_pair
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        real(DP), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: box_size(num_dimensions), max_over_box
        real(DP) :: permittivity

        if (interact) then
            domain%min = min_distance%get()
            data_field = prefix//"max distance over box edge"
            call input_data%get(data_field, max_over_box, data_found)
            call check_data_found(data_field, data_found)
            box_size = environment%periodic_box%get_size()
            domain%max = max_over_box * box_size(1)
            permittivity = environment%permittivity%get()
            select type (real_pair)
                type is (Tabulated_Ewald_Real_Pair)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call real_pair%construct(domain, permittivity, alpha)
    end subroutine construct_real_pair

    subroutine destroy_real_pair(real_pair)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(inout) :: real_pair

        if (allocated(real_pair)) then
            call real_pair%destroy()
            deallocate(real_pair)
        end if
    end subroutine destroy_real_pair

    subroutine allocate_and_construct_weighted_structure(weighted_structure)
        class(Abstract_Ewald_Reci_Structures), allocatable, intent(out) :: weighted_structure
    end subroutine allocate_and_construct_weighted_structure

    subroutine set_alpha(alpha, periodic_box, dipoles_exist, input_data, prefix)
        real(DP), intent(out) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: dipoles_exist
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: box_size(num_dimensions), alpha_times_box

        if (dipoles_exist) then
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

    subroutine set_are_dipolar(are_dipolar, components)
        logical, intent(out) :: are_dipolar(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(are_dipolar)
            are_dipolar(i_component) = component_is_dipolar(components(i_component)%dipolar_moments)
        end do
    end subroutine set_are_dipolar

    subroutine check_consistency(dipoles_exist, reciprocal_lattice, permittivity)
        logical, intent(in) :: dipoles_exist
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Permittivity), intent(in) :: permittivity

        if (dipoles_exist) then
            if (.not.use_reciprocal_lattice(reciprocal_lattice)) then
                call warning_continue("check_consistency: dipoles exist but reciprocal_lattice "//&
                    "unused.")
            end if
            if (.not.use_permittivity(permittivity)) then
                call warning_continue("check_consistency: dipoles exist but permittivity unused.")
            end if
        end if
    end subroutine

end module procedures_long_interactions_factory
