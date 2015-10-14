module procedures_metropolis_factory

use data_constants, only: num_components
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler, Concrete_Tower_Sampler, Null_Tower_Sampler
use class_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use types_metropolis_wrapper, only: Metropolis_Wrapper
use procedures_property_inquirers, only: component_can_move, component_can_rotate

implicit none

private
public :: metropolis_factory_create, metropolis_factory_set, metropolis_factory_destroy

interface metropolis_factory_create
    module procedure :: metropolis_factory_create_all
end interface metropolis_factory_create

interface metropolis_factory_set
    module procedure :: metropolis_factory_set_all
end interface metropolis_factory_set

interface metropolis_factory_destroy
    module procedure :: metropolis_factory_destroy_all
end interface metropolis_factory_destroy

contains

    subroutine metropolis_factory_create_all(metropolis, environment, changes)
        type(Metropolis_Wrapper), intent(out) :: metropolis
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        call allocate_and_construct_one_particle_move(metropolis%one_particle_move, environment, &
            changes)
        call allocate_and_construct_one_particle_rotation(metropolis%one_particle_rotation, &
            environment, changes)
    end subroutine metropolis_factory_create_all

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        class(Abstract_Tower_Sampler), allocatable :: selector

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

    subroutine metropolis_factory_set_all(metropolis, components, short_potentials, ewalds)
        type(Metropolis_Wrapper), intent(inout) :: metropolis
        type(Component_Wrapper), intent(in) :: components(num_components)
        type(Mixture_Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds

        call metropolis%one_particle_move%set_candidates(components, short_potentials, ewalds)
        call metropolis%one_particle_rotation%set_candidates(components, short_potentials, ewalds)
    end subroutine metropolis_factory_set_all

    subroutine allocate_and_construct_one_particle_rotation(one_particle_rotation, environment, &
        changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

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

    subroutine metropolis_factory_destroy_all(metropolis)
        type(Metropolis_Wrapper), intent(inout) :: metropolis

        call destroy_and_deallocate_one_particle_change(metropolis%one_particle_rotation)
        call destroy_and_deallocate_one_particle_change(metropolis%one_particle_move)
    end subroutine metropolis_factory_destroy_all

    subroutine destroy_and_deallocate_one_particle_change(one_particle_change)
        class(Abstract_One_Particle_Change), allocatable, intent(inout) :: one_particle_change

        call one_particle_change%destroy()
        if (allocated(one_particle_change)) deallocate(one_particle_change)
    end subroutine destroy_and_deallocate_one_particle_change

end module procedures_metropolis_factory
