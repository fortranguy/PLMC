module procedures_des_real_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use json_module, only: json_file
use procedures_checks, only: check_data_found
use class_periodic_box, only: Abstract_Periodic_Box
use class_permittivity, only: Abstract_Permittivity
use types_environment_wrapper, only: Environment_Wrapper
use class_minimum_distance, only: Abstract_Minimum_Distance
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Minimum_Distances_Wrapper
use types_potential_domain, only: Dipolar_Potential_Domain
use class_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use class_des_real_pair, only: Abstract_DES_Real_Pair, Tabulated_DES_Real_Pair, &
    Raw_DES_Real_Pair, Null_DES_Real_Pair
use class_des_real_component, only: Abstract_DES_Real_Component, &
    Concrete_DES_Real_Component, Null_DES_Real_Component
use class_des_real_visitor, only: Abstract_DES_Real_Visitor, Concrete_DES_Real_Visitor, &
    Null_DES_Real_Visitor
use types_dipolar_interactions_wrapper, only: DES_Real_Pairs_Wrapper, DES_Real_Component_Wrapper
use procedures_property_inquirers, only: components_interact

implicit none

private
public :: des_real_create, des_real_destroy

interface des_real_create
    module procedure :: create_visitor
    module procedure :: create_components
    module procedure :: create_component
    module procedure :: create_pairs
    module procedure :: create_pair
end interface des_real_create

interface des_real_destroy
    module procedure :: destroy_pair
    module procedure :: destroy_pairs
    module procedure :: destroy_component
    module procedure :: destroy_components
    module procedure :: destroy_visitor
end interface des_real_destroy

contains

    subroutine create_visitor(visitor, periodic_box, dipoles_exist)
        class(Abstract_DES_Real_Visitor), allocatable, intent(out) :: visitor
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            allocate(Concrete_DES_Real_Visitor :: visitor)
        else
            allocate(Null_DES_Real_Visitor :: visitor)
        end if
        call visitor%construct(periodic_box)
    end subroutine create_visitor

    subroutine destroy_visitor(visitor)
        class(Abstract_DES_Real_Visitor), allocatable, intent(inout) :: visitor

        if (allocated(visitor)) then
            call visitor%destroy()
            deallocate(visitor)
        end if
    end subroutine destroy_visitor

    subroutine create_components(components, periodic_box, mixture_components, real_pairs)
        type(DES_Real_Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        type(DES_Real_Pairs_Wrapper), intent(in) :: real_pairs(:)

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        allocate(components(size(real_pairs), size(real_pairs)))

        do j_component = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                associate (potential_ij => real_pairs(j_pair)%line(i_pair)%potential)
                    call des_real_create(components(i_component, j_component)%&
                        component, periodic_box, mixture_components(i_component)%positions, &
                        mixture_components(i_component)%dipolar_moments, potential_ij)
                end associate
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(DES_Real_Component_Wrapper), allocatable, intent(inout) :: components(:, :)

        integer :: i_component, j_component

        if (allocated(components)) then
            do j_component = size(components, 2), 1, -1
                do i_component = size(components, 1), 1, -1
                    call des_real_destroy(components(i_component, j_component)%component)
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine create_component(component, periodic_box, positions, dipolar_moments, pair)
        class(Abstract_DES_Real_Component), allocatable, intent(out) :: component
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        class(Abstract_DES_Real_Pair), intent(in) :: pair

        if (components_interact(pair)) then
            allocate(Concrete_DES_Real_Component :: component)
        else
            allocate(Null_DES_Real_Component :: component)
        end if
        call component%construct(periodic_box, positions, dipolar_moments, pair)
    end subroutine create_component

    subroutine destroy_component(component)
        class(Abstract_DES_Real_Component), allocatable, intent(inout) :: component

        if (allocated(component)) then
            call component%destroy()
            deallocate(component)
        end if
    end subroutine destroy_component

    subroutine create_pairs(real_pairs, permittivity, min_distances, are_dipolar, alpha, &
        input_data, prefix)
        type(DES_Real_Pairs_Wrapper), allocatable, intent(out) :: real_pairs(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Minimum_Distances_Wrapper), intent(in) :: min_distances(:)
        logical, intent(in) :: are_dipolar(:)
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        logical :: interact_ij
        integer :: i_component, j_component

        allocate(real_pairs(size(are_dipolar)))
        do j_component = 1, size(are_dipolar)
            allocate(real_pairs(j_component)%line(j_component))
            do i_component = 1, size(real_pairs(j_component)%line)
                interact_ij = are_dipolar(i_component) .and. are_dipolar(j_component)
                associate (min_distance_ij => min_distances(j_component)%line(i_component)%distance)
                    call des_real_create(real_pairs(j_component)%line(i_component)%potential, &
                        permittivity, min_distance_ij, interact_ij, alpha, input_data, prefix)
                end associate
            end do
        end do
    end subroutine create_pairs

    subroutine destroy_pairs(real_pairs)
        type(DES_Real_Pairs_Wrapper), allocatable, intent(inout) :: real_pairs(:)

        integer :: i_component, j_component

        if (allocated(real_pairs)) then
            do j_component = size(real_pairs), 1, -1
                if (allocated(real_pairs(j_component)%line)) then
                    do i_component = size(real_pairs(j_component)%line), 1, -1
                        call des_real_destroy(real_pairs(j_component)%line(i_component)%potential)
                    end do
                    deallocate(real_pairs(j_component)%line)
                end if
            end do
            deallocate(real_pairs)
        end if
    end subroutine destroy_pairs

    subroutine create_pair(pair, permittivity, min_distance, interact, alpha, input_data, prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call allocate_pair(pair, interact, input_data, prefix)
        call construct_pair(pair, permittivity, min_distance, interact, alpha, input_data, prefix)
    end subroutine create_pair

    subroutine allocate_pair(pair, interact, input_data, prefix)
        class(Abstract_DES_Real_Pair), allocatable, intent(out) :: pair
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
                allocate(Tabulated_DES_Real_Pair :: pair)
            else
                allocate(Raw_DES_Real_Pair :: pair)
            end if
        else
            allocate(Null_DES_Real_Pair :: pair)
        end if
    end subroutine allocate_pair

    subroutine construct_pair(pair, permittivity, min_distance, interact, alpha,  input_data, &
        prefix)
        class(Abstract_DES_Real_Pair), intent(inout) :: pair
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Minimum_Distance), intent(in) :: min_distance
        logical, intent(in) :: interact
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        type(Dipolar_Potential_Domain) :: domain

        if (interact) then
            domain%min = min_distance%get()
            data_field = prefix//"max distance over box edge"
            call input_data%get(data_field, domain%max_over_box, data_found)
            call check_data_found(data_field, data_found)
            select type (pair)
                type is (Tabulated_DES_Real_Pair)
                    data_field = prefix//"delta distance"
                    call input_data%get(data_field, domain%delta, data_found)
                    call check_data_found(data_field, data_found)
            end select
        end if
        call pair%construct(permittivity, alpha, domain)
    end subroutine construct_pair

    subroutine destroy_pair(pair)
        class(Abstract_DES_Real_Pair), allocatable, intent(inout) :: pair

        if (allocated(pair)) then
            call pair%destroy()
            deallocate(pair)
        end if
    end subroutine destroy_pair

end module procedures_des_real_factory
