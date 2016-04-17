module procedures_metropolis_algorithms_factory

use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use class_hetero_couples, only: Abstract_Hetero_Couples, Null_Hetero_Couples, &
    Concrete_Hetero_Couples
use types_component_wrapper, only: Component_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use class_two_particles_switch, only: Abstract_Two_Particles_Switch, &
    Concrete_Two_Particles_Switch, Null_Two_Particles_Switch
use types_metropolis_algorithms_wrapper, only: Metropolis_Algorithms_Wrapper
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: metropolis_algorithms_create, metropolis_algorithms_set, metropolis_algorithms_destroy, &
    metropolis_algorithms_create_move, metropolis_algorithms_create_rotation, &
    metropolis_algorithms_set_move, metropolis_algorithms_set_rotation


interface metropolis_algorithms_create
    module procedure :: create_all
    module procedure :: create_switch
end interface metropolis_algorithms_create

interface metropolis_algorithms_set
    module procedure :: set_all
    module procedure :: set_switch
end interface metropolis_algorithms_set

interface metropolis_algorithms_destroy
    module procedure :: destroy_switch
    module procedure :: destroy_change
    module procedure :: destroy_all
end interface metropolis_algorithms_destroy

contains

    subroutine create_all(metropolis_algorithms, environment, mixture, changes, short_interactions,&
        dipolar_interactions)
        type(Metropolis_Algorithms_Wrapper), intent(out) :: metropolis_algorithms
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        call metropolis_algorithms_create_move(metropolis_algorithms%one_particle_move, &
            environment, mixture, changes, short_interactions, dipolar_interactions)
        call metropolis_algorithms_create_rotation(metropolis_algorithms%one_particle_rotation, &
            environment, mixture, changes, short_interactions, dipolar_interactions)
        call metropolis_algorithms_create(metropolis_algorithms%two_particles_switch, environment, &
            mixture%components, short_interactions, dipolar_interactions)
    end subroutine create_all

    subroutine destroy_all(metropolis_algorithms)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms

        call metropolis_algorithms_destroy(metropolis_algorithms%two_particles_switch)
        call metropolis_algorithms_destroy(metropolis_algorithms%one_particle_rotation)
        call metropolis_algorithms_destroy(metropolis_algorithms%one_particle_move)
    end subroutine destroy_all

    subroutine metropolis_algorithms_create_move(one_particle_move, environment, mixture, &
        changes, short_interactions, dipolar_interactions)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        logical :: some_components_can_move
        integer :: i_component

        some_components_can_move = .false.
        do i_component = 1, size(changes)
            if (component_can_move(changes(i_component)%moved_positions)) then
                some_components_can_move = .true.
                exit
            end if
        end do

        if (some_components_can_move) then
            allocate(Concrete_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_One_Particle_Change :: one_particle_move)
        end if

        call one_particle_move%construct(environment, mixture, changes, short_interactions, &
            dipolar_interactions)
    end subroutine metropolis_algorithms_create_move

    subroutine metropolis_algorithms_create_rotation(one_particle_rotation, environment, mixture, &
        changes, short_interactions, dipolar_interactions)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Mixture_Wrapper), intent(in) :: mixture
        type(Changes_Component_Wrapper), intent(in) :: changes(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        logical :: some_components_can_rotate
        integer :: i_component

        some_components_can_rotate = .false.
        do i_component = 1, size(changes)
            if (component_can_rotate(changes(i_component)%rotated_orientations)) then
                some_components_can_rotate = .true.
                exit
            end if
        end do

        if (some_components_can_rotate) then
            allocate(Concrete_One_Particle_Rotation :: one_particle_rotation)
        else
            allocate(Null_One_Particle_Change :: one_particle_rotation)
        end if

        call one_particle_rotation%construct(environment, mixture, changes, short_interactions, &
            dipolar_interactions)
    end subroutine metropolis_algorithms_create_rotation

    subroutine destroy_change(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        if (allocated(one_particle_change)) then
            call one_particle_change%destroy()
            deallocate(one_particle_change)
        end if
    end subroutine destroy_change

    subroutine create_switch(two_particles_switch, environment, components, short_interactions, &
        dipolar_interactions)
        class(Abstract_Two_Particles_Switch), allocatable, intent(out) :: two_particles_switch
        type(Environment_Wrapper), intent(in) :: environment
        type(Component_Wrapper), intent(in) :: components(:)
        type(Short_Interactions_Wrapper), intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), intent(in) :: dipolar_interactions

        if (size(components) > 1) then
            allocate(Concrete_Two_Particles_Switch :: two_particles_switch)
        else
            allocate(Null_Two_Particles_Switch :: two_particles_switch)
        end if

        call two_particles_switch%construct(environment, components, short_interactions, &
            dipolar_interactions)
    end subroutine create_switch

    subroutine destroy_switch(two_particles_switch)
        class(Abstract_Two_Particles_Switch), allocatable, intent(inout) :: two_particles_switch

        if (allocated(two_particles_switch)) then
            call two_particles_switch%destroy()
            deallocate(two_particles_switch)
        end if
    end subroutine destroy_switch

    subroutine set_all(metropolis_algorithms, components, changes)
        type(Metropolis_Algorithms_Wrapper), intent(inout) :: metropolis_algorithms
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        call metropolis_algorithms_set_move(metropolis_algorithms%one_particle_move, components, &
            changes)
        call metropolis_algorithms_set_rotation(metropolis_algorithms%one_particle_rotation, &
            components, changes)
        call metropolis_algorithms_set(metropolis_algorithms%two_particles_switch, components)
    end subroutine set_all

    subroutine metropolis_algorithms_set_move(one_particle_move, components, changes)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_move
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        integer :: nums_candidates(size(components)), i_component

        do i_component = 1, size(nums_candidates)
            if (component_can_move(changes(i_component)%moved_positions)) then
                nums_candidates(i_component) = components(i_component)%average_number%get()
            else
                nums_candidates(i_component) = 0
            end if
        end do
        call allocate_selector(one_particle_move, nums_candidates)
    end subroutine metropolis_algorithms_set_move

    subroutine metropolis_algorithms_set_rotation(one_particle_rotation, components, changes)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_rotation
        type(Component_Wrapper), intent(in) :: components(:)
        type(Changes_Component_Wrapper), intent(in) :: changes(:)

        integer :: nums_candidates(size(components)), i_component

        do i_component = 1, size(nums_candidates)
            if (component_can_rotate(changes(i_component)%rotated_orientations)) then
                nums_candidates(i_component) = components(i_component)%average_number%get()
            else
                nums_candidates(i_component) = 0
            end if
        end do
        call allocate_selector(one_particle_rotation, nums_candidates)
    end subroutine metropolis_algorithms_set_rotation

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

    subroutine set_switch(two_particles_switch, components)
        class(Abstract_Two_Particles_Switch), intent(inout) :: two_particles_switch
        type(Component_Wrapper), intent(in) :: components(:)

        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector

        integer, allocatable :: nums_candidates(:)
        integer :: i_candidate, ij_couple(2)

        if (size(components) > 1) then
            allocate(Concrete_Hetero_Couples :: couples)
            allocate(Concrete_Tower_Sampler :: selector)
        else
            allocate(Null_Hetero_Couples :: couples)
            allocate(Null_Tower_Sampler :: selector)
        end if
        call couples%construct(size(components))

        allocate(nums_candidates(couples%get_num_indices()))
        do i_candidate = 1, size(nums_candidates)
            ij_couple = couples%get(i_candidate)
            nums_candidates(i_candidate) = minval([components(ij_couple(1))%number%get(), &
                components(ij_couple(2))%number%get()])
        end do
        call selector%construct(nums_candidates)

        call two_particles_switch%allocate_couples_and_selector(couples, selector)
        call selector%destroy()
        call couples%destroy()
    end subroutine set_switch

end module procedures_metropolis_algorithms_factory
