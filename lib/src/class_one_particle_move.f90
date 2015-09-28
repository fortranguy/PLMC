module class_one_particle_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_temperature, only: Abstract_Temperature
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_moved_positions, only: Abstract_Moved_Positions
use class_visitable_cells, only: Abstract_Visitable_Cells
use module_particle_energy, only: Concrete_Particle_Energy, &
    particle_energy_sum => Concrete_Particle_Energy_sum, operator(-)
use procedures_random, only: random_integer

implicit none

private

    type, abstract, public :: Abstract_One_Particle_Move
    private
        class(Abstract_Temperature), pointer :: temperature
        class(Abstract_Particles_Positions), pointer :: actor_positions, spectator_positions
        class(Abstract_Moved_Positions), pointer :: actor_moved_positions
        class(Abstract_Visitable_Cells), pointer :: actor_visitable_cells, inter_visitable_cells
    contains
        procedure :: construct => Abstract_One_Particle_Move_construct
        procedure :: destroy => Abstract_One_Particle_Move_destroy
        procedure :: set_actor => Abstract_One_Particle_Move_set_actor
        procedure :: set_spectator => Abstract_One_Particle_Move_set_spectator
        procedure :: try => Abstract_One_Particle_Move_try
    end type Abstract_One_Particle_Move

contains

    subroutine Abstract_One_Particle_Move_construct(this, temperature)
        class(Abstract_One_Particle_Move), intent(out) :: this
        class(Abstract_Temperature), target, intent(in) :: temperature

        this%temperature => temperature
    end subroutine Abstract_One_Particle_Move_construct

    subroutine Abstract_One_Particle_Move_destroy(this)
        class(Abstract_One_Particle_Move), intent(inout) :: this

        this%temperature => null()
    end subroutine Abstract_One_Particle_Move_destroy

    subroutine Abstract_One_Particle_Move_set_actor(this, actor_positions, actor_moved_positions, &
            actor_visitable_cells)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        class(Abstract_Particles_Positions), target, intent(in) :: actor_positions
        class(Abstract_Moved_Positions), target, intent(in) :: actor_moved_positions
        class(Abstract_Visitable_Cells), target, intent(in) :: actor_visitable_cells

        this%actor_positions => actor_positions
        this%actor_moved_positions => actor_moved_positions
        this%actor_visitable_cells => actor_visitable_cells
    end subroutine Abstract_One_Particle_Move_set_actor

    subroutine Abstract_One_Particle_Move_set_spectator(this, spectator_positions, &
        inter_visitable_cells)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        class(Abstract_Particles_Positions), target, intent(in) :: spectator_positions
        class(Abstract_Visitable_Cells), target, intent(in) :: inter_visitable_cells

        this%spectator_positions => spectator_positions
        this%inter_visitable_cells => inter_visitable_cells
    end subroutine Abstract_One_Particle_Move_set_spectator

    subroutine Abstract_One_Particle_Move_try(this, success, energy_difference)
        class(Abstract_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Particle_Energy), intent(out) :: energy_difference

        type(Concrete_Particle) :: old, new
        type(Concrete_Particle_Energy) :: new_energy, old_energy
        real(DP) :: energy_difference_sum
        logical :: overlap
        real(DP) :: rand

        old%i = random_integer(this%actor_positions%get_num())
        old%position = this%actor_positions%get(old%i)
        new%i = old%i
        new%position = this%actor_moved_positions%get(new%i)
        success = .false.
        if (this%actor_positions%get_num() > this%spectator_positions%get_num()) then
            call this%actor_visitable_cells%visit(overlap, new_energy%intra, new)
            if (overlap) return !Where?
            call this%inter_visitable_cells%visit(overlap, new_energy%inter, new)
        else
            call this%inter_visitable_cells%visit(overlap, new_energy%inter, new)
            if (overlap) return !Where?
            call this%actor_visitable_cells%visit(overlap, new_energy%intra, new)
        end if
        if (overlap) return !Where?
        call this%actor_visitable_cells%visit(overlap, old_energy%intra, old)
        call this%inter_visitable_cells%visit(overlap, old_energy%inter, old)

        energy_difference = new_energy - old_energy
        energy_difference_sum = particle_energy_sum(energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference_sum/this%temperature%get())) then
            call this%actor_positions%set(new%i, new%position)
            call this%actor_visitable_cells%move(old, new)
            call this%inter_visitable_cells%move(old, new)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Move_try

end module class_one_particle_move
