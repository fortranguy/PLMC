module procedures_metropolis_algorithms_factory

use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: metropolis_algorithms_create, metropolis_algorithms_set, metropolis_algorithms_destroy, &
    metropolis_algorithms_create_move, metropolis_algorithms_create_rotation


interface metropolis_algorithms_create
    module procedure :: create_all
end interface metropolis_algorithms_create

interface metropolis_algorithms_set
    module procedure :: set_all
end interface metropolis_algorithms_set

interface metropolis_algorithms_destroy
    module procedure :: destroy_change
    module procedure :: destroy_all
end interface metropolis_algorithms_destroy

contains

    subroutine create_all(metropolis_algorithms, environment, changes)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis_algorithms
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        call metropolis_algorithms_create_move(metropolis_algorithms%one_particle_move, &
            environment, changes)
        call metropolis_algorithms_create_rotation(metropolis_algorithms%one_particle_rotation, &
            environment, changes)
    end subroutine create_all

    subroutine destroy_all(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call metropolis_algorithms_destroy(metropolis_algorithms%one_particle_rotation)
        call metropolis_algorithms_destroy(metropolis_algorithms%one_particle_move)
    end subroutine destroy_all

    subroutine metropolis_algorithms_create_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: nums_candidates(size(changes))
        logical :: some_components_can_move
        integer :: i_component

        some_components_can_move = .false.
        do i_component = 1, size(nums_candidates)
            some_components_can_move = some_components_can_move .or. &
                component_can_move(changes(i_component)%moved_positions)
            nums_candidates(i_component) = changes(i_component)%moved_positions%get_num()
        end do
        if (some_components_can_move) then
            allocate(Concrete_Tower_Sampler :: selector)
            allocate(Concrete_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_Tower_Sampler :: selector)
            allocate(Null_One_Particle_Change :: one_particle_move)
        end if
        call selector%construct(nums_candidates)
        call one_particle_move%construct(environment, changes, selector)
        call selector%destroy()
    end subroutine metropolis_algorithms_create_move

    subroutine metropolis_algorithms_create_rotation(one_particle_rotation, environment, &
        changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        class(Abstract_Tower_Sampler), allocatable :: selector
        integer :: nums_candidates(size(changes))
        logical :: some_components_can_rotate
        integer :: i_component

        some_components_can_rotate = .false.
        do i_component = 1, size(changes)
            some_components_can_rotate = some_components_can_rotate .or. &
                component_can_rotate(changes(i_component)%rotated_orientations)
            nums_candidates(i_component) = changes(i_component)%rotated_orientations%get_num()
        end do
        if (some_components_can_rotate) then
            allocate(Concrete_Tower_Sampler :: selector)
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_Tower_Sampler :: selector)
            allocate(Null_One_Particle_Change :: one_particle_rotation)
        end if
        call selector%construct(nums_candidates)
        call one_particle_rotation%construct(environment, changes, selector)
        call selector%destroy()
    end subroutine metropolis_algorithms_create_rotation

    subroutine destroy_change(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        if (allocated(one_particle_change)) then
            call one_particle_change%destroy()
            deallocate(one_particle_change)
        end if
    end subroutine destroy_change

    subroutine set_all(metropolis_algorithms, components, short_interactions, long_interactions)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), intent(in) :: long_interactions

        call metropolis_algorithms%one_particle_move%set_candidates(components, &
            short_interactions, long_interactions)
        call metropolis_algorithms%one_particle_rotation%set_candidates(components, &
            short_interactions, long_interactions)
    end subroutine set_all

end module procedures_metropolis_algorithms_factory
