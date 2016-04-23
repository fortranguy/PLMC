module procedures_des_real_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use types_component_wrapper, only: Component_Wrapper
use types_min_distance_wrapper, only: Min_Distances_Line
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter
use types_des_real_pair_wrapper, only: DES_Real_Pairs_Line
use procedures_des_real_pair_factory, only: des_real_pair_create => create, &
    des_real_pair_destroy => destroy
use types_des_real_component_wrapper, only: DES_Real_Component_Wrapper
use procedures_des_real_component_factory, only: des_real_component_create => create, &
    des_real_component_destroy => destroy
use procedures_des_real_visitor_factory, only: des_real_visitor_create => create, &
    des_real_visitor_destroy => destroy

implicit none

private
public :: create, destroy

interface create
    module procedure :: des_real_visitor_create
    module procedure :: create_components
    module procedure :: des_real_component_create
    module procedure :: create_pairs
    module procedure :: des_real_pair_create
end interface create

interface destroy
    module procedure :: des_real_pair_destroy
    module procedure :: destroy_pairs
    module procedure :: des_real_component_destroy
    module procedure :: destroy_components
    module procedure :: des_real_visitor_destroy
end interface destroy

contains

    subroutine create_components(components, periodic_box, mixture_components, real_pairs)
        type(DES_Real_Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(Component_Wrapper), intent(in) :: mixture_components(:)
        type(DES_Real_Pairs_Line), intent(in) :: real_pairs(:)

        integer :: i_component, j_component
        integer :: i_pair, j_pair

        allocate(components(size(real_pairs), size(real_pairs)))

        do j_component = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                j_pair = maxval([i_component, j_component])
                i_pair = minval([i_component, j_component])
                associate (potential_ij => real_pairs(j_pair)%line(i_pair)%potential)
                    call create(components(i_component, j_component)%component, &
                        periodic_box, mixture_components(i_component)%positions, &
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
                    call destroy(components(i_component, j_component)%component)
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine create_pairs(real_pairs, permittivity, min_distances, are_dipolar, alpha, &
        input_data, prefix)
        type(DES_Real_Pairs_Line), allocatable, intent(out) :: real_pairs(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Min_Distances_Line), intent(in) :: min_distances(:)
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
                    call create(real_pairs(j_component)%line(i_component)%potential, &
                        permittivity, min_distance_ij, interact_ij, alpha, input_data, prefix)
                end associate
            end do
        end do
    end subroutine create_pairs

    subroutine destroy_pairs(real_pairs)
        type(DES_Real_Pairs_Line), allocatable, intent(inout) :: real_pairs(:)

        integer :: i_component, j_component

        if (allocated(real_pairs)) then
            do j_component = size(real_pairs), 1, -1
                if (allocated(real_pairs(j_component)%line)) then
                    do i_component = size(real_pairs(j_component)%line), 1, -1
                        call destroy(real_pairs(j_component)%line(i_component)%potential)
                    end do
                    deallocate(real_pairs(j_component)%line)
                end if
            end do
            deallocate(real_pairs)
        end if
    end subroutine destroy_pairs

end module procedures_des_real_factory
