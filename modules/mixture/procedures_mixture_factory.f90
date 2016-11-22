module procedures_mixture_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_input_prefixes, only: mixture_prefix
use json_module, only: json_file
use classes_number_to_string, only: Concrete_Number_to_String
use procedures_checks, only: check_data_found, check_positive
use classes_periodic_box, only: Abstract_Periodic_Box
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use procedures_component_factory, only: component_create => create, component_destroy => destroy
use procedures_composition_factory, only: composition_create => create, composition_destroy => &
    destroy
use procedures_hard_core_factory, only: hard_core_create => create, hard_core_destroy => destroy
use procedures_mixture_total_moments_factory, only: mixture_total_moments_create => create, &
    mixture_total_moments_destroy => destroy
use types_mixture_wrapper, only: Mixture_Wrapper
use procedures_mixture_inquirers, only: property_num_components => num_components, &
    mixture_can_exchange

implicit none

private
public :: create, destroy, rescale_positions

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

        logical :: components_exist, can_exchange

        call create(mixture%components, components_exist, can_exchange, environment%periodic_boxes,&
            generating_data, mixture_prefix)
        call composition_create(mixture%average_nums_particles, environment%&
            accessible_domains, mixture%components, components_exist, can_exchange)
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
        call composition_destroy(mixture%average_nums_particles)
        call destroy(mixture%components)
    end subroutine destroy_all

    subroutine create_components(components, components_exist, can_exchange, periodic_boxes, &
        generating_data, prefix)
        type(Component_Wrapper), allocatable, intent(out) :: components(:, :)
        logical, intent(out) :: components_exist, can_exchange
        class(Abstract_Periodic_Box), intent(in) :: periodic_boxes(:)
        type(json_file), intent(inout) :: generating_data
        character(len=*), intent(in) :: prefix

        integer :: i_box
        integer :: num_components, i_component
        type(Concrete_Number_to_String) :: string

        num_components = property_num_components(generating_data, prefix)
        call check_positive("create_components", "num_components", num_components)
        if (num_components == 0) then
            components_exist = .false.
            num_components = 1 !null component
            can_exchange = .false.
        else
            components_exist = .true.
            can_exchange = mixture_can_exchange(generating_data, prefix)
        end if
        allocate(components(num_components, size(periodic_boxes)))

        do i_box = 1, size(components, 2)
            do i_component = 1, size(components, 1)
                call component_create(components(i_component, i_box), periodic_boxes(i_box), &
                    components_exist, can_exchange, generating_data, prefix//"Component "//string%&
                        get(i_component)//".")
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

    subroutine rescale_positions(components, box_size_ratio)
        type(Component_Wrapper), intent(inout) :: components(:)
        real(DP), intent(in) :: box_size_ratio(:)

        integer :: i_component

        do i_component = 1, size(components)
            call components(i_component)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine rescale_positions

end module procedures_mixture_factory
