module procedures_metropolis_factory

use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use types_metropolis_wrapper, only: Metropolis_Wrapper
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: metropolis_create, metropolis_set, metropolis_destroy

interface metropolis_create
    module procedure :: metropolis_create_all
end interface metropolis_create

interface metropolis_set
    module procedure :: metropolis_set_all
end interface metropolis_set

interface metropolis_destroy
    module procedure :: metropolis_destroy_all
end interface metropolis_destroy

contains

    subroutine metropolis_create_all(metropolis, environment, changes)
        type(Metropolis_Wrapper), intent(out) :: metropolis
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(:)

        call allocate_and_construct_one_particle_move(metropolis%one_particle_move, environment, &
            changes)
        call allocate_and_construct_one_particle_rotation(metropolis%one_particle_rotation, &
            environment, changes)
    end subroutine metropolis_create_all

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        ! to update: multi components
        if (component_can_move(changes(1)%moved_positions) .or. &
            component_can_move(changes(2)%moved_positions)) then
            allocate(Concrete_Tower_Sampler :: selector)
            allocate(Concrete_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_Tower_Sampler :: selector)
            allocate(Null_One_Particle_Change :: one_particle_move)
        end if
        call one_particle_move%construct(environment, changes, selector)
        deallocate(selector)
    end subroutine allocate_and_construct_one_particle_move

    subroutine metropolis_set_all(metropolis, components, short_interactions, long_interactions)
        type(Metropolis_Wrapper), intent(inout) :: metropolis
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        call metropolis%one_particle_move%set_candidates(components, short_interactions, &
            long_interactions)
        call metropolis%one_particle_rotation%set_candidates(components, short_interactions, &
            long_interactions)
    end subroutine metropolis_set_all

    subroutine allocate_and_construct_one_particle_rotation(one_particle_rotation, environment, &
        changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(:)

        class(Abstract_Tower_Sampler), allocatable :: selector

        if (component_can_rotate(changes(1)%rotated_orientations) .or. &
            component_can_rotate(changes(2)%rotated_orientations)) then
            allocate(Concrete_Tower_Sampler :: selector)
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_Tower_Sampler :: selector)
            allocate(Null_One_Particle_Change :: one_particle_rotation)
        end if
        call one_particle_rotation%construct(environment, changes, selector)
        deallocate(selector)
    end subroutine allocate_and_construct_one_particle_rotation

    subroutine metropolis_destroy_all(metropolis)
        type(Metropolis_Wrapper), intent(inout) :: metropolis

        call destroy_and_deallocate_one_particle_change(metropolis%one_particle_rotation)
        call destroy_and_deallocate_one_particle_change(metropolis%one_particle_move)
    end subroutine metropolis_destroy_all

    subroutine destroy_and_deallocate_one_particle_change(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        call one_particle_change%destroy()
        if (allocated(one_particle_change)) deallocate(one_particle_change)
    end subroutine destroy_and_deallocate_one_particle_change

end module procedures_metropolis_factory
