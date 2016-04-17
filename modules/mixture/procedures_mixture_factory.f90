module procedures_mixture_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use procedures_checks, only: check_data_found, check_positive
use class_number_to_string, only: Concrete_Number_to_String
use class_periodic_box, only: Abstract_Periodic_Box
use class_walls_potential, only: Abstract_Walls_Potential
use types_environment_wrapper, only: Environment_Wrapper
use class_mixture_total_moment, only: Abstract_Mixture_Total_Moment, &
    Concrete_Mixture_Total_Moment, Null_Mixture_Total_Moment
use types_component_wrapper, only: Component_Wrapper
use class_minimum_distance, only: Abstract_Minimum_Distance, Concrete_Minimum_Distance, &
    Null_Minimum_Distance
use procedures_component_factory, only: component_create, component_destroy
use types_mixture_wrapper, only: Minimum_Distance_Wrapper, Minimum_Distances_Wrapper, &
    Mixture_Wrapper
use types_readers_wrapper, only: Component_Readers_wrapper
use procedures_property_inquirers, only: use_walls, component_exists, component_is_dipolar

implicit none

private
public :: mixture_create, mixture_destroy, mixture_set

interface mixture_create
    module procedure :: create_all
    module procedure :: create_components
    module procedure :: create_components_min_distances
    module procedure :: create_min_distance
    module procedure :: create_wall_min_distances
    module procedure :: create_total_moment
end interface mixture_create

interface mixture_destroy
    module procedure :: destroy_total_moment
    module procedure :: destroy_min_distance
    module procedure :: destroy_min_distances
    module procedure :: destroy_components_min_distances
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface mixture_destroy

interface mixture_set
    module procedure :: set_are_dipolar
    module procedure :: read_mixture_coordinates
end interface mixture_set

contains

    subroutine create_all(mixture, environment, input_data, prefix)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(Environment_Wrapper), intent(in) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call mixture_create(mixture%components, environment%periodic_box, input_data, prefix)
        call mixture_create(mixture%components_min_distances, mixture%components, input_data, &
            prefix)
        call mixture_create(mixture%wall_min_distances, mixture%components, environment%&
            walls_potential, input_data, prefix)
        call mixture_create(mixture%total_moment, mixture%components)
    end subroutine create_all

    subroutine destroy_all(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        call mixture_destroy(mixture%total_moment)
        call mixture_destroy(mixture%wall_min_distances)
        call mixture_destroy(mixture%components_min_distances)
        call mixture_destroy(mixture%components)
    end subroutine destroy_all

    subroutine create_components(components, periodic_box, input_data, prefix)
        type(Component_Wrapper), allocatable, intent(out) :: components(:)
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: exists, data_found
        integer :: num_components, i_component
        type(Concrete_Number_to_String) :: string

        data_field = prefix//"number of components"
        call input_data%get(data_field, num_components, data_found)
        call check_data_found(data_field, data_found)
        call check_positive("create_components", "num_components", num_components)
        if (num_components == 0) then
            exists = .false.
            num_components = 1 !null component
        else
            exists = .true.
        end if
        allocate(components(num_components))
        do i_component = 1, size(components)
            call component_create(components(i_component), exists, periodic_box, &
                input_data, prefix//"Component "//string%get(i_component)//".")
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(Component_Wrapper), allocatable, intent(inout) :: components(:)

        integer :: i_component

        if (allocated(components)) then
            do i_component = size(components), 1, -1
                call component_destroy(components(i_component))
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine create_components_min_distances(min_distances, components, input_data, prefix)
        type(Minimum_Distances_Wrapper), allocatable, intent(out) :: min_distances(:)
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
                call mixture_create(min_distances(j_component)%line(i_component)%distance, exists, &
                    input_data, min_distance_prefix)
                deallocate(min_distance_prefix)
            end do
        end do
    end subroutine create_components_min_distances

    subroutine destroy_components_min_distances(min_distances)
        type(Minimum_Distances_Wrapper), allocatable, intent(inout) :: min_distances(:)

        integer :: i_component

        if (allocated(min_distances)) then
            do i_component = size(min_distances), 1, -1
                call mixture_destroy(min_distances(i_component)%line)
            end do
            deallocate(min_distances)
        end if
    end subroutine destroy_components_min_distances

    subroutine create_wall_min_distances(min_distances, components, potential, input_data, prefix)
        type(Minimum_Distance_Wrapper), allocatable, intent(out) :: min_distances(:)
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
            call mixture_create(min_distances(i_component)%distance, exists, input_data, &
                prefix//"Component "//string%get(i_component)//".With Walls.")
        end do
    end subroutine create_wall_min_distances

    subroutine destroy_min_distances(min_distances)
        type(Minimum_Distance_Wrapper), allocatable, intent(inout) :: min_distances(:)

        integer :: i_component

        if (allocated(min_distances)) then
            do i_component = size(min_distances), 1, -1
                call mixture_destroy(min_distances(i_component)%distance)
            end do
            deallocate(min_distances)
        end if
    end subroutine destroy_min_distances

    subroutine create_min_distance(min_distance, exists, input_data, prefix)
        class(Abstract_Minimum_Distance), allocatable, intent(out) :: min_distance
        logical, intent(in) :: exists
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        character(len=:), allocatable :: data_field
        logical :: data_found
        real(DP) :: min_distance_value

        if (exists) then
            allocate(Concrete_Minimum_Distance :: min_distance)
            data_field = prefix//"minimum distance"
            call input_data%get(data_field, min_distance_value, data_found)
            call check_data_found(data_field, data_found)
        else
            allocate(Null_Minimum_Distance :: min_distance)
            min_distance_value = 0._DP
        end if
        call min_distance%set(min_distance_value)
    end subroutine create_min_distance

    subroutine destroy_min_distance(min_distance)
        class(Abstract_Minimum_Distance), allocatable, intent(inout) :: min_distance

        if (allocated(min_distance)) deallocate(min_distance)
    end subroutine destroy_min_distance

    subroutine create_total_moment(total_moment, components)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(out) :: total_moment
        type(Component_Wrapper), intent(in) :: components(:)

        logical :: are_dipolar(size(components))

        call mixture_set(are_dipolar, components)
        if (any(are_dipolar)) then
            allocate(Concrete_Mixture_Total_Moment :: total_moment)
        else
            allocate(Null_Mixture_Total_Moment :: total_moment)
        end if
        call total_moment%construct(components, are_dipolar)
    end subroutine create_total_moment

    subroutine destroy_total_moment(total_moment)
        class(Abstract_Mixture_Total_Moment), allocatable, intent(inout) :: total_moment

        if (allocated(total_moment)) then
            call total_moment%destroy()
            deallocate(total_moment)
        end if
    end subroutine destroy_total_moment

    subroutine set_are_dipolar(are_dipolar, components)
        logical, intent(out) :: are_dipolar(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(are_dipolar)
            are_dipolar(i_component) = component_is_dipolar(components(i_component)%dipolar_moments)
        end do
    end subroutine set_are_dipolar

    subroutine read_mixture_coordinates(components_readers, input_data, prefix)
        type(Component_Readers_wrapper), intent(in) :: components_readers(:)
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        integer :: i_component
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field, filename
        logical :: data_found

        do i_component = 1, size(components_readers)
            data_field = prefix//"Component "//string%get(i_component)//"."//"initial coordinates"
            call input_data%get(data_field, filename, data_found)
            call check_data_found(data_field, data_found)
            call components_readers(i_component)%coordinates%read(filename)
        end do
    end subroutine read_mixture_coordinates

end module procedures_mixture_factory
