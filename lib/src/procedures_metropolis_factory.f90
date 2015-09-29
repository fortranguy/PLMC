module procedures_metropolis_factory

use json_module, only: json_file
use module_data, only: test_data_found
use class_temperature, only: Abstract_Temperature
use class_moved_positions, only: Abstract_Moved_Positions
use class_one_particle_move, only: Abstract_One_Particle_Move, &
    Null_One_Particle_Move, Two_Candidates_One_Particle_Move, &
    First_Candidate_One_Particle_Move, Second_Candidate_One_Particle_Move
use procedures_types_selectors, only: particles_can_move

implicit none

private
public :: metropolis_factory_create, metropolis_factory_destroy

interface metropolis_factory_create
    module procedure :: allocate_and_construct_one_particle_move
end interface metropolis_factory_create

interface metropolis_factory_destroy
    module procedure :: destroy_and_deallocate_one_particle_move
end interface metropolis_factory_destroy

contains

    subroutine allocate_and_construct_one_particle_move(one_particle_move, temperature, &
        moved_positions_1, moved_positions_2)
        class(Abstract_One_Particle_Move), allocatable, intent(out) :: one_particle_move
        class(Abstract_Temperature), intent(in) :: temperature
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
        call one_particle_move%construct(temperature, moved_positions_1, moved_positions_2)
    end subroutine allocate_and_construct_one_particle_move

    subroutine destroy_and_deallocate_one_particle_move(one_particle_move)
        class(Abstract_One_Particle_Move), allocatable, intent(inout) :: one_particle_move

        call one_particle_move%destroy()
        if (allocated(one_particle_move)) deallocate(one_particle_move)
    end subroutine destroy_and_deallocate_one_particle_move

end module procedures_metropolis_factory
