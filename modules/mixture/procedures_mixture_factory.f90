module procedures_mixture_factory

use data_input_prefixes, only: mixture_prefix
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found, check_positive
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_component_factory, only: component_create => create, component_destroy => destroy
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use procedures_mixture_total_moments_factory, only: mixture_total_moments_create => create, &
    mixture_total_moments_destroy => destroy
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_inquirers, only: component_has_positions, component_has_orientations, &
    property_num_components => num_components, mixture_can_exchange

implicit none

private
public :: create, destroy, set_nums_particles, set_have_positions, &
    set_have_orientations

interface create
    module procedure :: create_all
    module procedure :: create_components
end interface create

interface destroy
    module procedure :: destroy_components
    module procedure :: destroy_all
end interface destroy

contains

    !> @todo more rigorous components(:) input for components_min_distances and wall_min_distances?
    subroutine create_all(mixture, environment, generating_data)
        type(Mixture_Wrapper), intent(out) :: mixture
        type(Environment_Wrapper), intent(in) :: environment
        type(json_file), intent(inout) :: generating_data

        call create(mixture%components, environment%periodic_boxes, environment%&
            accessible_domains, generating_data, mixture_prefix)
        call hard_core_create(mixture%components_min_distances, mixture%components(:, 1), &
            generating_data, mixture_prefix)
        call hard_core_create(mixture%wall_min_distances, environment%wall_min_distance, mixture%&
            components_min_distances, mixture%components(:, 1), environment%visitable_walls(1))
        call mixture_total_moments_create(mixture%total_moments, mixture%components)
    end subroutine create_all

    subroutine destroy_all(mixture)
        type(Mixture_Wrapper), intent(inout) :: mixture

        call mixture_total_moments_destroy(mixture%total_moments)
        call hard_core_destroy(mixture%wall_min_distances)
        call hard_core_destroy(mixture%components_min_distances)
        call destroy(mixture%components)
    end subroutine destroy_all

    subroutine create_components(components, periodic_boxes, accessible_domains, generating_data, &
        prefix)
        type(Component_Wrapper), allocatable, intent(out) :: components(:, :)
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        class(Abstract_Parallelepiped_Domain), intent(in) :: accessible_domains(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        logical :: exists, can_exchange
        integer :: num_components, i_component
        type(Concrete_Number_to_String) :: string

        num_components = property_num_components(generating_data, prefix)
        call check_positive("create_components", "num_components", num_components)
        if (num_components == 0) then
            exists = .false.
            num_components = 1 !null component
            can_exchange = .false.
        else
            exists = .true.
            can_exchange = mixture_can_exchange(generating_data, prefix)
        end if
        allocate(components(num_components, size(periodic_boxes)))

        do i_box = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                call component_create(components(i_component, i_box), periodic_boxes(i_box), &
                    size(periodic_boxes) == 1, accessible_domains(i_box), exists, can_exchange, &
                    generating_data, prefix//"Component "//string%get(i_component)//".")
            end do
        end do
    end subroutine create_components

    subroutine destroy_components(components)
        type(Component_Wrapper), allocatable, intent(inout) :: components(:, :)

        integer :: i_box, i_component

        if (allocated(components)) then
            do i_box = size(components, 2), 1, -1
                do i_component = size(components, 1), 1, -1
                    call component_destroy(components(i_component, i_box))
                end do
            end do
            deallocate(components)
        end if
    end subroutine destroy_components

    subroutine set_nums_particles(nums_particles, components)
        integer, intent(inout) :: nums_particles(:)
        type(Component_Wrapper), intent(in) :: components(:)

        integer :: i_component

        do i_component = 1, size(nums_particles)
            nums_particles(i_component) = components(i_component)%num_particles%get()
        end do
    end subroutine set_nums_particles

    subroutine set_have_positions(have_positions, components)
        logical, intent(inout) :: have_positions(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(have_positions, 2)
            do i_component = 1, size(have_positions, 1)
                have_positions(i_component, i_box) = &
                    component_has_positions(components(i_component, i_box)%positions)
            end do
        end do
    end subroutine set_have_positions

    subroutine set_have_orientations(have_orientations, components)
        logical, intent(inout) :: have_orientations(:, :)
        type(Component_Wrapper), intent(in) :: components(:, :)

        integer :: i_box, i_component

        do i_box = 1, size(have_orientations, 2)
            do i_component = 1, size(have_orientations, 1)
                have_orientations(i_component, i_box) = &
                    component_has_orientations(components(i_component, i_box)%orientations)
            end do
        end do
    end subroutine set_have_orientations

end module procedures_mixture_factory
