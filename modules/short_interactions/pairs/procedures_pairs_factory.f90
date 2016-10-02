module procedures_pairs_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use types_min_distance_wrapper, only: Min_Distance_Wrapper, Min_Distances_Line
use classes_potential_expression, only: Abstract_Potential_Expression
use procedures_potential_expression_factory, only: potential_expression_create => create, &
    potential_expression_destroy => destroy
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use procedures_pair_potential_factory, only: pair_potential_create => create, &
    pair_potential_destroy => destroy
use procedures_short_pairs_visitors_factory, only: short_pairs_visitors_create => create, &
    short_pairs_visitors_destroy => destroy
use procedures_short_interactions_inquirers, only: components_interact

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_wall
    module procedure :: create_components
    module procedure :: short_pairs_visitors_create
    module procedure :: pair_potential_create
    module procedure :: potential_expression_create
end interface create

interface destroy
    module procedure :: potential_expression_destroy
    module procedure :: pair_potential_destroy
    module procedure :: short_pairs_visitors_destroy
    module procedure :: destroy_components
    module procedure :: destroy_pairs
end interface

contains

    subroutine create_wall(pairs, interact, min_distances, generating_data, prefix)
        type(Pair_Potential_Wrapper), allocatable, intent(out) :: pairs(:)
        logical, intent(out) :: interact
        type(Min_Distance_Wrapper), intent(in) :: min_distances(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        class(Abstract_Potential_Expression), allocatable :: expression
        character(len=:), allocatable :: pair_prefix
        type(Concrete_Number_to_String) :: string

        interact = .true.
        allocate(pairs(size(min_distances)))
        do i_component = 1, size(pairs)
            pair_prefix = prefix//"Component "//string%get(i_component)//".With Walls."
            associate (min_distance => min_distances(i_component)%distance, &
                interact_i => components_interact(min_distances(i_component)%distance))
                interact = interact .and. interact_i
                call create(expression, interact_i, generating_data, pair_prefix)
                call create(pairs(i_component)%potential, min_distance, expression, interact_i, &
                    generating_data, pair_prefix)
            end associate
            call destroy(expression)
        end do
    end subroutine create_wall

    subroutine create_components(pairs, interact, min_distances, generating_data, prefix)
        type(Pair_Potentials_Line), allocatable, intent(out) :: pairs(:)
        logical, intent(out) :: interact
        type(Min_Distances_Line), intent(in) :: min_distances(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_component, j_component
        logical :: interact_ij
        class(Abstract_Potential_Expression), allocatable :: expression
        character(len=:), allocatable :: pair_prefix
        type(Concrete_Number_to_String) :: string

        interact = .true.
        allocate(pairs(size(min_distances)))
        do j_component = 1, size(pairs)
            allocate(pairs(j_component)%line(j_component))
            do i_component = 1, size(pairs(j_component)%line)
                if (i_component == j_component) then
                    pair_prefix = prefix//"Component "//string%get(i_component)//"."
                else
                    pair_prefix = prefix//"Inter "//string%get(i_component)//&
                        string%get(j_component)//"."
                end if
                associate (min_distance => min_distances(j_component)%line(i_component)%distance)
                    interact_ij = components_interact(min_distance)
                    interact = interact .and. interact_ij
                    call create(expression, interact_ij, generating_data, pair_prefix)
                    call create(pairs(j_component)%line(i_component)%potential, min_distance, &
                        expression, interact_ij, generating_data, pair_prefix)
                end associate
                call destroy(expression)
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(pairs)
        type(Pair_Potentials_Line), allocatable, intent(inout) :: pairs(:)

        integer :: i_component

        if (allocated(pairs)) then
            do i_component = size(pairs), 1, -1
                call destroy(pairs(i_component)%line)
            end do
            deallocate(pairs)
        end if
    end subroutine destroy_components

    subroutine destroy_pairs(pairs)
        type(Pair_Potential_Wrapper), allocatable, intent(inout) :: pairs(:)

        integer :: i_component

        if (allocated(pairs)) then
            do i_component = size(pairs), 1, -1
                call destroy(pairs(i_component)%potential)
            end do
            deallocate(pairs)
        end if
    end subroutine destroy_pairs

end module procedures_pairs_factory
