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
use types_component_wrapper, only: Component_Wrapper, Mixture_Wrapper_Old
use types_mixture_wrapper, only: Minimum_Distances_Wrapper, Are_Dipolar_Wrapper
use types_potential_domain, only: Concrete_Potential_Domain
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair, &
    Tabulated_Ewald_Real_Pair, Raw_Ewald_Real_Pair, Null_Ewald_Real_Pair
use class_ewald_real_component, only: Abstract_Ewald_Real_Component, &
    Concrete_Ewald_Real_Component, Null_Ewald_Real_Component
use class_weighted_structure, only: Abstract_Weighted_Structure
use types_ewalds_wrapper, only: Ewald_Wrapper, Ewald_Wrapper_Macro, Ewald_Wrapper_Micro, &
    Ewald_Real_Pair_Wrapper, Ewald_Real_Pairs_Wrapper, Ewalds_Wrapper
use procedures_property_inquirers, only: components_interact

implicit none

private
public :: ewald_create, ewald_destroy

interface ewald_create
    module procedure :: ewald_create_all
    module procedure :: create_real_pairs
    module procedure :: ewald_create_macro
    module procedure :: ewald_create_micro
    module procedure :: allocate_and_construct_real_pair
    module procedure :: allocate_and_construct_real_component
end interface ewald_create

interface ewald_destroy
    module procedure :: destroy_and_deallocate_real_component
    module procedure :: destroy_and_deallocate_real_pair
    module procedure :: destroy_real_pairs
    module procedure :: ewald_destroy_micro
    module procedure :: ewald_destroy_macro
    module procedure :: ewald_destroy_all
end interface ewald_destroy

contains

    subroutine ewald_create_all(ewald, environment, component, input_data, prefix)
        type(Ewald_Wrapper), intent(out) :: ewald
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: component
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: alpha
        logical :: are_dipolar

        !are_dipolar = component_are_dipolar(component%dipolar_moments)
        call set_alpha(alpha, are_dipolar, environment%periodic_box, input_data, prefix)
        !call ewald_create(ewald%real_pair, are_dipolar, alpha, environment%periodic_box, &
        !    component%diameter, input_data, prefix//"Real.")
        call ewald_create(ewald%real_component, are_dipolar, environment%periodic_box, &
            component%positions, component%dipolar_moments, ewald%real_pair)
    end subroutine ewald_create_all

    subroutine ewald_destroy_all(ewald)
        type(Ewald_Wrapper), intent(inout) :: ewald

        call ewald_destroy(ewald%real_component)
        call ewald_destroy(ewald%real_pair)
    end subroutine ewald_destroy_all

    subroutine create_real_pairs(real_pairs, periodic_box, min_distances, are_dipolar, alpha, &
        input_data, prefix)
        type(Ewald_Real_Pairs_Wrapper), allocatable, intent(out) :: real_pairs(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Minimum_Distances_Wrapper), intent(in) :: min_distances(:)
        type(Are_Dipolar_Wrapper), intent(in) :: are_dipolar(:)
        real(DP), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: j_component, i_component

        allocate(real_pairs(size(are_dipolar)))
        do j_component = 1, size(are_dipolar)
            allocate(real_pairs(j_component)%with_components(j_component))
            do i_component = 1, size(real_pairs(j_component)%with_components)
                associate (are_dipolar_ij => are_dipolar(j_component)%&
                    with_components(i_component), &
                    min_distance => min_distances(j_component)%with_components(i_component)%&
                    min_distance)
                    call ewald_create(real_pairs(j_component)%with_components(i_component)%&
                        real_pair, are_dipolar_ij, alpha, periodic_box, min_distance, input_data, &
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

    subroutine ewald_create_macro(ewald_macro, ewald_micro, environment, component)
        type(Ewald_Wrapper_Macro), intent(out) :: ewald_macro
        type(Ewald_Wrapper_Micro), intent(in) :: ewald_micro
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: component

        logical :: interact

        interact = components_interact(ewald_micro%real_pair)
        call ewald_create(ewald_macro%real_component, interact, environment%periodic_box, &
            component%positions, component%dipolar_moments, ewald_micro%real_pair)
    end subroutine ewald_create_macro

    subroutine ewald_create_micro(ewald_micro, are_dipolar, environment, min_distance, &
        input_data, prefix)
        type(Ewald_Wrapper_Micro), intent(out) :: ewald_micro
        logical, intent(in) :: are_dipolar
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        real(DP) :: alpha

        call set_alpha(alpha, are_dipolar, environment%periodic_box, input_data, prefix)
        call ewald_create(ewald_micro%real_pair, are_dipolar, alpha, &
            environment%periodic_box, min_distance, input_data, prefix//"Real.")
    end subroutine ewald_create_micro

    subroutine set_alpha(alpha, are_dipolar, periodic_box, input_data, prefix)
        real(DP), intent(out) :: alpha
        logical, intent(in) :: are_dipolar
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: box_size(num_dimensions), alpha_times_box

        if (are_dipolar) then
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

    subroutine allocate_and_construct_real_pair(real_pair, are_dipolar, alpha, periodic_box, &
        min_distance, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: are_dipolar
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_real_pair(real_pair, are_dipolar, input_data, prefix)
        call construct_real_pair(real_pair, are_dipolar, alpha, periodic_box, min_distance, &
            input_data, prefix)
    end subroutine allocate_and_construct_real_pair

    subroutine allocate_real_pair(real_pair, are_dipolar, input_data, prefix)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(out) :: real_pair
        logical, intent(in) :: are_dipolar
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found, tabulated_potential

        if (are_dipolar) then
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

    subroutine construct_real_pair(real_pair, are_dipolar, alpha, periodic_box, min_distance, &
        input_data, prefix)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: real_pair
        logical, intent(in) :: are_dipolar
        real(DP), intent(in) :: alpha
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: box_size(num_dimensions), max_over_box

        if (are_dipolar) then
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

    subroutine allocate_and_construct_real_component(real_component, are_dipolar, periodic_box, &
        component_positions, component_dipolar_moments, real_pair)
        class(Abstract_Ewald_Real_Component), allocatable, intent(out) :: real_component
        logical, intent(in) :: are_dipolar
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: component_dipolar_moments
        class(Abstract_Ewald_Real_Pair), intent(in) :: real_pair

        if (are_dipolar) then
            allocate(Concrete_Ewald_Real_Component :: real_component)
        else
            allocate(Null_Ewald_Real_Component :: real_component)
        end if
        call real_component%construct(periodic_box, component_positions, &
            component_dipolar_moments, real_pair)
    end subroutine allocate_and_construct_real_component

    subroutine allocate_and_construct_weighted_structure(weighted_structure)
        class(Abstract_Weighted_Structure), allocatable, intent(out) :: weighted_structure
    end subroutine allocate_and_construct_weighted_structure

    subroutine ewald_destroy_micro(ewald_micro)
        type(Ewald_Wrapper_Micro), intent(inout) :: ewald_micro

        call ewald_destroy(ewald_micro%real_pair)
    end subroutine ewald_destroy_micro

    subroutine ewald_destroy_macro(ewald_macro)
        type(Ewald_Wrapper_Macro), intent(inout) :: ewald_macro

        call ewald_destroy(ewald_macro%real_component)
    end subroutine ewald_destroy_macro

    subroutine destroy_and_deallocate_real_component(real_component)
        class(Abstract_Ewald_Real_Component), allocatable, intent(inout) :: real_component

        call real_component%destroy()
        if (allocated(real_component)) deallocate(real_component)
    end subroutine destroy_and_deallocate_real_component

    subroutine destroy_and_deallocate_real_pair(real_pair)
        class(Abstract_Ewald_Real_Pair), allocatable, intent(inout) :: real_pair

        call real_pair%destroy()
        if (allocated(real_pair)) deallocate(real_pair)
    end subroutine destroy_and_deallocate_real_pair

end module procedures_ewalds_factory
