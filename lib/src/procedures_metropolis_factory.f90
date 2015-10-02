module procedures_metropolis_factory

use json_module, only: json_file
use module_data, only: test_data_found
use procedures_property_inquirers, only: particles_can_move
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Particles_Wrapper
use types_changes_wrapper, only: Changes_Wrapper
use types_short_potential_wrapper, only: Short_Potential_Wrapper, Short_Potential_Macro_Wrapper, &
    Wall_Short_Potential_Wrapper
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
    module procedure :: set_one_particle_move_components
    module procedure :: set_one_particle_move_observables
    module procedure :: set_one_particle_move_short_potentials
end interface metropolis_factory_set

interface metropolis_factory_destroy
    module procedure :: destroy_and_deallocate_one_particle_move
end interface metropolis_factory_destroy

contains

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, changes)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        type(Changes_Wrapper), intent(in) :: changes(2)

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

    subroutine set_one_particle_move_components(one_particle_move, components)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Particles_Wrapper), intent(in) :: components(2)

        call one_particle_move%set_candidate(1, components(1)%positions)
        call one_particle_move%set_candidate(2, components(2)%positions)
    end subroutine set_one_particle_move_components

    subroutine set_one_particle_move_short_potentials(one_particle_move, intras, inters, walls)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Short_Potential_Wrapper), intent(in) :: intras(2)
        type(Short_Potential_Macro_Wrapper), intent(in) :: inters(2)
        type(Wall_Short_Potential_Wrapper), intent(in) :: walls(2)

        call one_particle_move%set_candidate(1, intras(1)%cells, inters(1)%cells)
        call one_particle_move%set_candidate(1, walls(1)%pair)
        call one_particle_move%set_candidate(2, intras(2)%cells, inters(2)%cells)
        call one_particle_move%set_candidate(2, walls(2)%pair)
    end subroutine set_one_particle_move_short_potentials

    subroutine set_one_particle_move_observables(one_particle_move, observables)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Concrete_Mixture_Observables), intent(in) :: observables

        call one_particle_move%set_candidates_observables(observables%move_counters, &
            observables%particles_energies, observables%inter_energy)
    end subroutine set_one_particle_move_observables

    subroutine destroy_and_deallocate_one_particle_move(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        call one_particle_move%destroy()
        if (allocated(one_particle_move)) deallocate(one_particle_move)
    end subroutine destroy_and_deallocate_one_particle_move

end module procedures_metropolis_factory
