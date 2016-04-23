module procedures_mixture_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found, check_positive
use classes_number_to_string, only: Concrete_Number_to_String
use classes_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_component_factory, only: component_create => create, component_destroy => destroy
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use procedures_mixture_total_moment_factory, only: mixture_total_moment_create => create, &
    mixture_total_moment_destroy => destroy
use types_mixture_wrapper, only: Mixture_Wrapper
use types_readers_wrapper, only: Component_Coordinates_Reader_wrapper

implicit none

private
public :: mixture_create, mixture_destroy, mixture_set_initial_coordinates

interface mixture_create
    module procedure :: create_all
    module procedure :: create_components
end interface mixture_create

interface mixture_destroy
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface mixture_destroy

contains

    subroutine create_all(mixture, environment, input_data, prefix)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(Environment_Wrapper), intent(in) :: environment
        type(json_file), intent(inout) :: input_data
        character(len=*), intent(in) :: prefix

        call mixture_create(mixture%components, environment%periodic_box, input_data, prefix)
        call hard_core_create(mixture%wall_min_distances, mixture%components, environment%&
            walls, input_data, prefix)
        call hard_core_create(mixture%components_min_distances, mixture%components, input_data, &
            prefix)
        call mixture_total_moment_create(mixture%total_moment, mixture%components)
    end subroutine create_all

    subroutine destroy_all(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        call mixture_total_moment_destroy(mixture%total_moment)
        call hard_core_destroy(mixture%components_min_distances)
        call hard_core_destroy(mixture%wall_min_distances)
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

    subroutine mixture_set_initial_coordinates(components_readers, input_data, prefix)
        type(Component_Coordinates_Reader_wrapper), intent(in) :: components_readers(:)
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
    end subroutine mixture_set_initial_coordinates

end module procedures_mixture_factory
