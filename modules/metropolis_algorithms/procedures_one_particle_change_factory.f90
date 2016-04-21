module procedures_one_particle_change_factory

use classes_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use classes_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: create_move, create_rotation, destroy, set_move, set_rotation

contains

    subroutine create_move(one_particle_move, environment, mixture, &
        change_components, short_interactions, dipolar_interactions)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        logical :: some_components_can_move
        integer :: i_component

        some_components_can_move = .false.
        do i_component = 1, size(change_components)
            if (component_can_move(change_components(i_component)%moved_positions)) then
                some_components_can_move = .true.
                exit
            end if
        end do

        if (some_components_can_move) then
            allocate(Concrete_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_One_Particle_Change :: one_particle_move)
        end if

        call one_particle_move%construct(environment, mixture, change_components, &
            short_interactions, dipolar_interactions)
    end subroutine create_move

    subroutine create_rotation(one_particle_rotation, environment, mixture, &
        change_components, short_interactions, dipolar_interactions)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        logical :: some_components_can_rotate
        integer :: i_component

        some_components_can_rotate = .false.
        do i_component = 1, size(change_components)
            if (component_can_rotate(change_components(i_component)%rotated_orientations)) then
                some_components_can_rotate = .true.
                exit
            end if
        end do

        if (some_components_can_rotate) then
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_One_Particle_Change :: one_particle_rotation)
        end if

        call one_particle_rotation%construct(environment, mixture, change_components, &
            short_interactions, dipolar_interactions)
    end subroutine create_rotation

    subroutine destroy(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        if (allocated(one_particle_change)) then
            call one_particle_change%destroy()
            deallocate(one_particle_change)
        end if
    end subroutine destroy

    subroutine set_move(one_particle_move, components, change_components)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_move
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        integer :: nums_candidates(size(components)), i_component

        do i_component = 1, size(nums_candidates)
            if (component_can_move(change_components(i_component)%moved_positions)) then
                nums_candidates(i_component) = components(i_component)%average_number%get()
            else
                nums_candidates(i_component) = 0
            end if
        end do
        call allocate_selector(one_particle_move, nums_candidates)
    end subroutine set_move

    subroutine set_rotation(one_particle_rotation, components, change_components)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_rotation
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: change_components(:)

        integer :: nums_candidates(size(components)), i_component

        do i_component = 1, size(nums_candidates)
            if (component_can_rotate(change_components(i_component)%rotated_orientations)) then
                nums_candidates(i_component) = components(i_component)%average_number%get()
            else
                nums_candidates(i_component) = 0
            end if
        end do
        call allocate_selector(one_particle_rotation, nums_candidates)
    end subroutine set_rotation

    subroutine allocate_selector(one_particle_change, nums_candidates)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_change
        integer, intent(in) :: nums_candidates(:)

        class(Abstract_Tower_Sampler), allocatable :: selector

        if (any(nums_candidates /= 0)) then
            allocate(Concrete_Tower_Sampler :: selector)
        else
            allocate(Null_Tower_Sampler :: selector)
        end if

        call selector%construct(nums_candidates)
        call one_particle_change%allocate_selector(selector)
        call selector%destroy()
    end subroutine allocate_selector

end module procedures_one_particle_change_factory
