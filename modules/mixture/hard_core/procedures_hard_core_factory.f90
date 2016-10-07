module procedures_hard_core_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use classes_visitable_walls, only: Abstract_Visitable_Walls
use procedures_environment_inquirers, only: use_walls
use classes_min_distance, only: Abstract_Min_Distance
use types_component_wrapper, only: Component_Wrapper
use types_min_distance_wrapper, only: Min_Distance_Wrapper, Min_Distances_Line
use procedures_min_distance_factory, only: min_distance_create_from_json => create_from_json, &
    min_distance_create_from_value => create_from_value, min_distance_destroy => destroy
use procedures_mixture_inquirers, only: component_exists

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_wall_component
    module procedure :: create_components
    module procedure :: min_distance_create_from_json, min_distance_create_from_value
end interface

interface destroy
    module procedure :: min_distance_destroy
    module procedure :: destroy_components_min_distances
    module procedure :: destroy_min_distances
end interface

contains

    subroutine create_wall_component(min_distances, wall_min_distance, components_min_distances, &
        components, visitable_walls)
        type(Min_Distance_Wrapper), allocatable, intent(out) :: min_distances(:)
        class(Abstract_Min_Distance), intent(in) :: wall_min_distance
        type(Min_Distances_Line), intent(in) :: components_min_distances(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Visitable_Walls), intent(in) :: visitable_walls

        real(DP) :: min_distance_i
        integer :: i_component
        logical :: exists

        allocate(min_distances(size(components)))
        do i_component = 1, size(min_distances)
            exists = use_walls(visitable_walls) .and. component_exists(components(i_component)%&
                num_particles)
            min_distance_i = (wall_min_distance%get() + components_min_distances(i_component)%&
            line(i_component)%distance%get()) / 2._DP
            call create(min_distances(i_component)%distance, exists, min_distance_i)
        end do
    end subroutine create_wall_component

    subroutine destroy_min_distances(min_distances)
        type(Min_Distance_Wrapper), allocatable, intent(inout) :: min_distances(:)

        integer :: i_component

        if (allocated(min_distances)) then
            do i_component = size(min_distances), 1, -1
                call destroy(min_distances(i_component)%distance)
            end do
            deallocate(min_distances)
        end if
    end subroutine destroy_min_distances

     subroutine create_components(min_distances, components, generating_data, prefix)
        type(Min_Distances_Line), allocatable, intent(out) :: min_distances(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_component, j_component
        logical :: exists
        character(len=:), allocatable :: min_distance_prefix
        type(Concrete_Number_to_String) :: string

        allocate(min_distances(size(components)))
        do j_component = 1, size(min_distances)
            allocate(min_distances(j_component)%line(j_component))
            do i_component = 1, size(min_distances(j_component)%line)
                exists = component_exists(components(j_component)%num_particles) .and. &
                    component_exists(components(i_component)%num_particles)
                if (i_component == j_component) then
                    min_distance_prefix = prefix//"Component "//string%get(i_component)//"."
                else
                    min_distance_prefix = prefix//"Inter "//string%get(i_component)//&
                        string%get(j_component)//"."
                end if
                call create(min_distances(j_component)%line(i_component)%distance, exists, &
                    generating_data, min_distance_prefix)
                deallocate(min_distance_prefix)
            end do
        end do
    end subroutine create_components

    subroutine destroy_components_min_distances(min_distances)
        type(Min_Distances_Line), allocatable, intent(inout) :: min_distances(:)

        integer :: i_component

        if (allocated(min_distances)) then
            do i_component = size(min_distances), 1, -1
                call destroy(min_distances(i_component)%line)
            end do
            deallocate(min_distances)
        end if
    end subroutine destroy_components_min_distances

end module procedures_hard_core_factory
