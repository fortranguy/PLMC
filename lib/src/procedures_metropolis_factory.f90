module procedures_metropolis_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use procedures_property_inquirers, only: particles_can_move
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Particles_Wrapper
use class_moved_positions, only: Abstract_Moved_Positions
use types_short_potential_wrapper, only: Short_Potential_Wrapper, Short_Potential_Macro_Wrapper
use class_one_particle_move, only: Abstract_One_Particle_Move, &
    Null_One_Particle_Move, Two_Candidates_One_Particle_Move, &
    First_Candidate_One_Particle_Move, Second_Candidate_One_Particle_Move
use types_change_counter, only: Concrete_Change_Counter
use module_particle_energy, only: Concrete_Particle_Energy

implicit none

private
public :: metropolis_factory_create, metropolis_factory_set, metropolis_factory_destroy

interface metropolis_factory_create
    module procedure :: allocate_and_construct_one_particle_move
end interface metropolis_factory_create

interface metropolis_factory_set
    module procedure :: set_one_particle_move_components_and_potential
    module procedure :: set_one_particle_move_observables
end interface metropolis_factory_set

interface metropolis_factory_destroy
    module procedure :: destroy_and_deallocate_one_particle_move
end interface metropolis_factory_destroy

contains

    subroutine allocate_and_construct_one_particle_move(one_particle_move, environment, &
        moved_positions_1, moved_positions_2)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_move
        type(Environment_Wrapper), intent(in) :: environment
        class(Abstract_Moved_Positions), intent(in) :: moved_positions_1, moved_positions_2

        if (particles_can_move(moved_positions_1) .and. particles_can_move(moved_positions_2)) then
            allocate(Two_Candidates_One_Particle_Move :: one_particle_move)
        else if (particles_can_move(moved_positions_1) .and. &
            .not.particles_can_move(moved_positions_2)) then
            allocate(First_Candidate_One_Particle_Move :: one_particle_move)
        else if (.not.particles_can_move(moved_positions_1) .and. &
            particles_can_move(moved_positions_2)) then
            allocate(Second_Candidate_One_Particle_Move :: one_particle_move)
        else
            allocate(Null_One_Particle_Move :: one_particle_move)
        end if
        call one_particle_move%construct(environment, moved_positions_1, moved_positions_2)
    end subroutine allocate_and_construct_one_particle_move

    subroutine set_one_particle_move_components_and_potential(one_particle_move, components, &
        intras, inters)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Particles_Wrapper), intent(in) :: components(2)
        type(Short_Potential_Wrapper), intent(in) :: intras(2)
        type(Short_Potential_Macro_Wrapper), intent(in) :: inters(2)

        call one_particle_move%set_candidate(1, components(1)%positions)
        call one_particle_move%set_candidate(1, intras(1)%cells, inters(1)%cells)
        call one_particle_move%set_candidate(2, components(2)%positions)
        call one_particle_move%set_candidate(2, intras(2)%cells, inters(2)%cells)
    end subroutine set_one_particle_move_components_and_potential

    subroutine set_one_particle_move_observables(one_particle_move, move_counters, &
        particles_energies, inter_energy)
        class(Abstract_One_Particle_Move), intent(inout) :: one_particle_move
        type(Concrete_Change_Counter), intent(in) :: move_counters(2)
        type(Concrete_Particle_Energy), intent(in) :: particles_energies(2)
        real(DP), intent(in) :: inter_energy

        call one_particle_move%set_candidates_observables(move_counters, particles_energies, &
            inter_energy)
    end subroutine set_one_particle_move_observables

    subroutine destroy_and_deallocate_one_particle_move(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        call one_particle_move%destroy()
        if (allocated(one_particle_move)) deallocate(one_particle_move)
    end subroutine destroy_and_deallocate_one_particle_move

end module procedures_metropolis_factory
