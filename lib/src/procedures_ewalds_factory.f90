module procedures_ewalds_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use class_minimum_distance, only: Abstract_Minimum_Distance
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Minimum_Distances_Wrapper, Mixture_Wrapper
use types_potential_domain, only: Concrete_Potential_Domain
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, &
    Tabulated_Ewald_Real_Pair, Raw_Ewald_Real_Pair, Null_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component, &
    Concrete_Ewald_Real_Component, Null_Ewald_Real_Component
use class_weighted_structure, only: Abstract_Weighted_Structure
use types_ewalds_wrapper, only: Ewald_Wrapper, Ewald_Wrapper_Macro, Ewald_Wrapper_Micro, &
    Ewald_Real_Pair_Wrapper, Ewald_Real_Pairs_Wrapper, Ewald_Real_Component_Wrapper, Ewalds_Wrapper
use procedures_property_inquirers, only: components_interact

implicit none

private
public :: ewald_create, ewald_destroy

interface ewald_create
    module procedure :: create_all
    module procedure :: create_real_components
    module procedure :: create_real_component
    module procedure :: create_real_pairs
    module procedure :: create_real_pair
    module procedure :: set_alpha
end interface ewald_create

interface ewald_destroy
    module procedure :: destroy_real_pair
    module procedure :: destroy_real_pairs
    module procedure :: destroy_real_component
    module procedure :: destroy_real_components
    module procedure :: destroy_all
end interface ewald_destroy

contains

    subroutine create_all(ewalds, environment, mixture, input_data, prefix)
        type(Ewalds_Wrapper), intent(out) :: ewalds
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: alpha

        call ewald_create(ewalds%real_components, mixture%components_are_dipolar, &
            environment%periodic_box, mixture%components)
        call ewald_create(alpha, any(mixture%components_are_dipolar), environment%periodic_box, &
            input_data, prefix)
        call ewald_create(ewalds%real_pairs, mixture%components_are_dipolar, &
            environment%periodic_box, mixture%inter_min_distances, alpha, input_data, &
            prefix//"Real.")
    end subroutine create_all

    subroutine destroy_all(ewalds)
        type(Ewalds_Wrapper), intent(inout) :: ewalds

        call ewald_destroy(ewalds%real_pairs)
        call ewald_destroy(ewalds%real_components)
    end subroutine destroy_all

    subroutine create_real_components(real_components, components_are_dipolar, periodic_box, &
        components)
        type(Ewald_Real_Component_Wrapper), allocatable, intent(out) :: real_components(:)
        logical, intent(in) :: components_are_dipolar(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        allocate(real_components(size(components_are_dipolar)))
        do i_component = 1, size(real_components)
            call ewald_create(real_components(i_component)%real_component, &
                components_are_dipolar(i_component), periodic_box, &
                components(i_component)%positions, &
                components(i_component)%dipolar_moments)
        end do
    end subroutine create_real_components

    subroutine destroy_real_components(real_components)
        type(Ewald_Real_Component_Wrapper), allocatable, intent(inout) :: real_components(:)

        if (allocated(real_components)) deallocate(real_components)
    end subroutine destroy_real_components

    subroutine create_real_component(real_component, is_dipolar, periodic_box, &
        component_positions, component_dipolar_moments)
        class(Abstract_Ewald_Real_Component), allocatable, intent(out) :: real_component
        logical, intent(in) :: is_dipolar
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments

        if (is_dipolar) then
            allocate(Concrete_Ewald_Real_Component :: real_component)
        else
            allocate(Null_Ewald_Real_Component :: real_component)
        end if
        call real_component%construct(periodic_box, component_positions, component_dipolar_moments)
    end subroutine create_real_component

    subroutine destroy_real_component(real_component)
        class(Abstract_Ewald_Real_Component), allocatable, intent(inout) :: real_component

        call real_component%destroy()
        if (allocated(real_component)) deallocate(real_component)
    end subroutine destroy_real_component

    subroutine create_real_pairs(real_pairs, components_are_dipolar, periodic_box, min_distances, &
        alpha, input_data, prefix)
        type(Ewald_Real_Pairs_Wrapper), allocatable, intent(out) :: real_pairs(:)
        logical, intent(in) :: components_are_dipolar(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Minimum_Distances_Wrapper), intent(in) :: min_distances(:)
        real(DP), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: j_component, i_component

        allocate(real_pairs(size(components_are_dipolar)))
        do j_component = 1, size(components_are_dipolar)
            allocate(real_pairs(j_component)%with_components(j_component))
            do i_component = 1, size(real_pairs(j_component)%with_components)
                associate (interact_ij => components_are_dipolar(i_component) .and. &
                    components_are_dipolar(j_component), &
                    min_distance => min_distances(j_component)%with_components(i_component)%&
                    min_distance)
                    call ewald_create(real_pairs(j_component)%with_components(i_component)%&
                        real_pair, interact_ij, alpha, periodic_box, min_distance, input_data, &
                        prefix)
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
                        call ewald_destroy(real_pairs(j_component)%with_components(i_component)%&
                            real_pair)
                    end do
                    deallocate(real_pairs(j_component)%with_components)
                end if
            end do
            deallocate(real_pairs)
        end if
    end subroutine destroy_real_pairs

    subroutine set_alpha(alpha, dipoles_exist, periodic_box, input_data, prefix)
        real(DP), intent(out) :: alpha
        logical, intent(in) :: dipoles_exist
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
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

    subroutine create_real_pair(real_pair, interact, alpha, periodic_box, &
        min_distance, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: interact
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_real_pair(real_pair, interact, input_data, prefix)
        call construct_real_pair(real_pair, interact, alpha, periodic_box, min_distance, &
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

    subroutine construct_real_pair(real_pair, interact, alpha, periodic_box, min_distance, &
        input_data, prefix)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: real_pair
        logical, intent(in) :: interact
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: box_size(num_dimensions), max_over_box

        if (interact) then
            domain%min = min_distance%get()
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

    subroutine destroy_real_pair(real_pair)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(inout) :: real_pair

        call real_pair%destroy()
        if (allocated(real_pair)) deallocate(real_pair)
    end subroutine destroy_real_pair

    subroutine allocate_and_construct_weighted_structure(weighted_structure)
        class(Abstract_Weighted_Structure), allocatable, intent(out) :: weighted_structure
    end subroutine allocate_and_construct_weighted_structure

end module procedures_ewalds_factory
