module procedures_min_distances_factory

use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use classes_walls_potential, only: Abstract_Walls_Potential
use types_component_wrapper, only: Component_Wrapper
use procedures_min_distance_factory, only: min_distance_create => create, &
    min_distance_destroy => destroy
use types_min_distances_wrapper, only: Min_Distance_Wrapper, Min_Distances_Wrapper
use procedures_property_inquirers, only: use_walls, component_exists

implicit none

private
public :: create, destroy

interface create
    module procedure :: create_wall_min_distances
    module procedure :: create_components_min_distances
    module procedure :: min_distance_create
end interface

interface destroy
    module procedure :: min_distance_destroy
    module procedure :: destroy_components_min_distances
    module procedure :: destroy_min_distances
end interface

contains

    subroutine create_wall_min_distances(min_distances, components, potential, input_data, prefix)
        type(Min_Distance_Wrapper), allocatable, intent(out) :: min_distances(:)
        type(Component_Wrapper), intent(in) :: components(:)
        class(Abstract_Walls_Potential), intent(in) :: potential
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        logical :: exists
        type(Concrete_Number_to_String) :: string

        allocate(min_distances(size(components)))
        do i_component = 1, size(min_distances)
            exists = use_walls(potential) .and. component_exists(components(i_component)%number)
            call create(min_distances(i_component)%distance, exists, input_data, &
                prefix//"Component "//string%get(i_component)//".With Walls.")
        end do
    end subroutine create_wall_min_distances

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

     subroutine create_components_min_distances(min_distances, components, input_data, prefix)
        type(Min_Distances_Wrapper), allocatable, intent(out) :: min_distances(:)
        type(Component_Wrapper), intent(in) :: components(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_component, j_component
        logical :: exists
        character(len=:), allocatable :: min_distance_prefix
        type(Concrete_Number_to_String) :: string

        allocate(min_distances(size(components)))
        do j_component = 1, size(min_distances)
            allocate(min_distances(j_component)%line(j_component))
            do i_component = 1, size(min_distances(j_component)%line)
                exists = component_exists(components(j_component)%number) .and. &
                    component_exists(components(i_component)%number)
                if (i_component == j_component) then
                    min_distance_prefix = prefix//"Component "//string%get(i_component)//"."
                else
                    min_distance_prefix = prefix//"Inter "//string%get(i_component)//&
                        string%get(j_component)//"."
                end if
                call create(min_distances(j_component)%line(i_component)%distance, exists, &
                    input_data, min_distance_prefix)
                deallocate(min_distance_prefix)
            end do
        end do
    end subroutine create_components_min_distances

    subroutine destroy_components_min_distances(min_distances)
        type(Min_Distances_Wrapper), allocatable, intent(inout) :: min_distances(:)

        integer :: i_component

        if (allocated(min_distances)) then
            do i_component = size(min_distances), 1, -1
                call destroy(min_distances(i_component)%line)
            end do
            deallocate(min_distances)
        end if
    end subroutine destroy_components_min_distances

end module procedures_min_distances_factory
