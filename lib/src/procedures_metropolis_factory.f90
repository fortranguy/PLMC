module procedures_metropolis_factory

use data_constants, only: num_components
use json_module, only: json_file
use procedures_property_inquirers, only: particles_can_move
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Particles_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_one_particle_move, only: Abstract_One_Particle_Move, &
    Null_One_Particle_Move, Two_Candidates_One_Particle_Move, &
    First_Candidate_One_Particle_Move, Second_Candidate_One_Particle_Move
use types_mixture_observables, only: Concrete_Mixture_Observables

implicit none

private
public :: metropolis_factory_create, metropolis_factory_set, metropolis_factory_destroy

interface metropolis_factory_create
    module procedure :: allocate_and_construct_one_particle_move
end interface metropolis_factory_create

interface metropolis_factory_set
    module procedure :: set_one_particle_move
end interface metropolis_factory_set

interface metropolis_factory_destroy
    module procedure :: destroy_and_deallocate_one_particle_move
end interface metropolis_factory_destroy

contains

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(num_components)

        if (particles_can_move(changes(1)%moved_positions) .and. &
            particles_can_move(changes(2)%moved_positions)) then
            allocate(Two_Candidates_One_Particle_Move :: one_particle_move)
        else if (particles_can_move(changes(1)%moved_positions) .and. &
            .not.particles_can_move(changes(2)%moved_positions)) then
            allocate(First_Candidate_One_Particle_Move :: one_particle_move)
        else if (.not.particles_can_move(changes(1)%moved_positions) .and. &
            particles_can_move(changes(2)%moved_positions)) then
            allocate(Second_Candidate_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_One_Particle_Move :: one_particle_move)
        end if
        call one_particle_move%construct(environment, changes(1)%moved_positions, &
            changes(2)%moved_positions)
    end subroutine allocate_and_construct_one_particle_move

    subroutine set_one_particle_move(one_particle_move, components, short_potentials, ewalds, &
        observables)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Particles_Wrapper), intent(in) :: components(num_components)
        type(Mixture_Short_Potentials_Wrapper), intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), intent(in) :: ewalds
        type(Concrete_Mixture_Observables), intent(in) :: observables

        call one_particle_move%set(components)
        call one_particle_move%set(short_potentials)
        call one_particle_move%set(ewalds)
        call one_particle_move%set(observables%changes_counters%moves, &
            observables%particles_energies, observables%inter_energy)
    end subroutine set_one_particle_move

    subroutine destroy_and_deallocate_one_particle_move(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        call one_particle_move%destroy()
        if (allocated(one_particle_move)) deallocate(one_particle_move)
    end subroutine destroy_and_deallocate_one_particle_move

end module procedures_metropolis_factory
