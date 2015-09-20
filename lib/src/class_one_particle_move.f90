module class_one_particle_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_temperature, only: Abstract_Temperature
use types_particles, only: Concrete_Particle, Concrete_Particles
use class_moved_particles_positions, only: Abstract_Moved_Particles_Positions
use class_visitable_cells, only: Abstract_Visitable_Cells
use module_particle_energy, only: Concrete_Particle_Energy, &
    particle_energy_sum => Concrete_Particle_Energy_sum, operator(-)
use procedures_random, only: random_integer

implicit none

private

    type, abstract, public :: Abstract_One_Particle_Move
    private
        class(Abstract_Temperature), pointer :: temperature
        type(Concrete_Particles), pointer :: actor_particles, spectator_particles
        class(Abstract_Moved_Particles_Positions), pointer :: moved_positions
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

    subroutine Abstract_One_Particle_Move_set_actor(this, particles, moved_positions, &
            visitable_cells)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Concrete_Particles), target, intent(in) :: particles
        class(Abstract_Moved_Particles_Positions), target, intent(in) :: moved_positions
        class(Abstract_Visitable_Cells), target, intent(in) :: visitable_cells

        this%actor_particles => particles
        this%moved_positions => moved_positions
        this%actor_visitable_cells => visitable_cells
    end subroutine Abstract_One_Particle_Move_set_actor

    subroutine Abstract_One_Particle_Move_set_spectator(this, particles, inter_visitable_cells)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Concrete_Particles), target, intent(in) :: particles
        class(Abstract_Visitable_Cells), target, intent(in) :: inter_visitable_cells

        this%spectator_particles => particles
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

        old%i = random_integer(this%actor_particles%number%get())
        old%position = this%actor_particles%positions%get(old%i)
        new%i = old%i
        new%position = this%moved_positions%get(new%i)
        success = .false.
        if (this%actor_particles%number%get() > this%spectator_particles%number%get()) then
            call this%actor_visitable_cells%visit(new, overlap, new_energy%intra)
            if (overlap) return !Where?
            call this%inter_visitable_cells%visit(new, overlap, new_energy%inter)
        else
            call this%inter_visitable_cells%visit(new, overlap, new_energy%inter)
            if (overlap) return !Where?
            call this%actor_visitable_cells%visit(new, overlap, new_energy%intra)
        end if
        call this%actor_visitable_cells%visit(old, overlap, old_energy%intra)
        call this%inter_visitable_cells%visit(old, overlap, old_energy%inter)

        energy_difference = new_energy - old_energy
        energy_difference_sum = particle_energy_sum(energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference_sum/this%temperature%get())) then
            call this%actor_particles%positions%set(new%i, new%position)
            call this%actor_visitable_cells%move(old, new)
            call this%inter_visitable_cells%move(old, new)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Move_try

end module class_one_particle_move
