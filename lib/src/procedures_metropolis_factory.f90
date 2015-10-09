module procedures_metropolis_factory

use data_constants, only: num_components
use json_module, only: json_file
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Particles_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_candidates_selector, only: Abstract_Candidates_Selector, Concrete_Candidates_Selector, &
    First_Candidate_Selector, Second_Candidate_Selector, Null_Candidates_Selector
use class_one_particle_change, only: Abstract_One_Particle_Change, &
    Concrete_One_Particle_Move, Concrete_One_Particle_Rotation, Null_One_Particle_Change
use class_one_particle_move, only: Abstract_One_Particle_Move, Null_One_Particle_Move, &
    Two_Candidates_Move, First_Candidate_Move, Second_Candidate_Move
use class_one_particle_rotation, only: Abstract_One_Particle_Rotation, Null_One_Particle_Rotation, &
    Two_Candidates_Rotation, First_Candidate_Rotation, Second_Candidate_Rotation
use types_mixture_observables, only: Concrete_Mixture_Observables
use procedures_property_inquirers, only: particles_can_move, particles_can_rotate

implicit none

private
public :: metropolis_factory_create, metropolis_factory_set, metropolis_factory_destroy

interface metropolis_factory_create
    module procedure :: allocate_and_construct_one_particle_move
    module procedure :: allocate_and_construct_one_particle_rotation
end interface metropolis_factory_create

interface metropolis_factory_set
    module procedure :: set_one_particle_change
end interface metropolis_factory_set

interface metropolis_factory_destroy
    module procedure :: destroy_and_deallocate_one_particle_move
    module procedure :: destroy_and_deallocate_one_particle_rotation
end interface metropolis_factory_destroy

contains

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Change), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        class(Abstract_Candidates_Selector), allocatable :: selector

        if (particles_can_move(changes(1)%moved_positions) .and. &
            particles_can_move(changes(2)%moved_positions)) then
            allocate(Concrete_Candidates_Selector :: selector)
        else if (particles_can_move(changes(1)%moved_positions) .and. &
            .not.particles_can_move(changes(2)%moved_positions)) then
            allocate(First_Candidate_Selector :: selector)
        else if (.not.particles_can_move(changes(1)%moved_positions) .and. &
            particles_can_move(changes(2)%moved_positions)) then
            allocate(Second_Candidate_Selector :: selector)
        else
            allocate(Null_Candidates_Selector :: selector)
        end if
        select type (selector)
            type is (Null_Candidates_Selector)
                allocate(Null_One_Particle_Change :: one_particle_move)
            class default
                allocate(Concrete_One_Particle_Move :: one_particle_move)
        end select
        call one_particle_move%construct(environment, selector, [changes(1)%moved_positions, &
            changes(2)%moved_positions])
        if (allocated(selector)) deallocate(selector)
    end subroutine allocate_and_construct_one_particle_move

    subroutine set_one_particle_move(one_particle_change, components, short_potentials, ewalds, &
        observables)
        class(Concrete_One_Particle_Move), intent(inout) :: one_particle_change

    end subroutine set_one_particle_move

    subroutine set_one_particle_change(one_particle_change, components, short_potentials, ewalds, &
        observables)
        class(Abstract_One_Particle_Change), intent(inout) :: one_particle_change
        type(Particles_Wrapper), intent(in) :: components(num_components)
        type(Mixture_Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Concrete_Mixture_Observables), intent(in) :: observables

        call one_particle_change%set(components)
        call one_particle_change%set(short_potentials)
        call one_particle_change%set(ewalds)
        call one_particle_change%set(observables%changes_counters, &
            observables%particles_energies, observables%inter_energy)
    end subroutine set_one_particle_change

    subroutine allocate_and_construct_one_particle_rotation(one_particle_rotation, environment, &
        changes)
        class(Abstract_One_Particle_Rotation), allocatable, intent(out) :: one_particle_rotation
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        if (particles_can_rotate(changes(1)%rotated_orientations) .and. &
            particles_can_rotate(changes(2)%rotated_orientations)) then
            allocate(Two_Candidates_Rotation :: one_particle_rotation)
        else if (particles_can_rotate(changes(1)%rotated_orientations) .and. &
            .not.particles_can_rotate(changes(2)%rotated_orientations)) then
            allocate(First_Candidate_Rotation :: one_particle_rotation)
        else if (.not.particles_can_rotate(changes(1)%rotated_orientations) .and. &
            particles_can_rotate(changes(2)%rotated_orientations)) then
            allocate(Second_Candidate_Rotation :: one_particle_rotation)
        else
            allocate(Null_One_Particle_Rotation :: one_particle_rotation)
        end if
        call one_particle_rotation%construct(environment, changes(1)%rotated_orientations, &
            changes(2)%rotated_orientations)
    end subroutine allocate_and_construct_one_particle_rotation

    subroutine set_one_particle_rotation(one_particle_rotation, components, ewalds, observables)
        class(Abstract_One_Particle_Rotation), intent(inout) :: one_particle_rotation
        type(Particles_Wrapper), intent(in) :: components(num_components)
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Concrete_Mixture_Observables), intent(in) :: observables

        call one_particle_rotation%set(components)
        call one_particle_rotation%set(ewalds)
        call one_particle_rotation%set(observables%changes_counters%rotations, &
            observables%particles_energies, observables%inter_energy)
    end subroutine set_one_particle_rotation

    subroutine destroy_and_deallocate_one_particle_move(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        call one_particle_move%destroy()
        if (allocated(one_particle_move)) deallocate(one_particle_move)
    end subroutine destroy_and_deallocate_one_particle_move

    subroutine destroy_and_deallocate_one_particle_rotation(one_particle_rotation)
        class(Abstract_One_Particle_Rotation), allocatable, intent(inout) :: one_particle_rotation

        call one_particle_rotation%destroy()
        if (allocated(one_particle_rotation)) deallocate(one_particle_rotation)
    end subroutine destroy_and_deallocate_one_particle_rotation

end module procedures_metropolis_factory
