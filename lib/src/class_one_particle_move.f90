module class_one_particle_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components
use procedures_errors, only: error_exit
use procedures_checks, only: check_in_range
use procedures_random, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use class_temperature, only: Abstract_Temperature
use class_external_field, only: Abstract_External_Field
use class_walls_potential, only: Abstract_Walls_Potential
use types_particle, only: Concrete_Particle
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use types_particles_wrapper, only: Particles_Wrapper
use class_moved_positions, only: Abstract_Moved_Positions
use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_cells, only: Abstract_Visitable_Cells
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_particles, only: Abstract_Ewald_Real_Particles
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use module_particles_energy, only: Concrete_Particles_Energy, Concrete_Inter_Energy, &
    particle_energy_sum => Concrete_Particles_Energy_sum, &
    inter_energy_sum => Concrete_Inter_Energy_sum, operator(+), operator(-)
use module_changes_success, only: Concrete_Change_Counters

implicit none

private

    type :: Move_Candidate
        class(Abstract_Particles_Positions), pointer :: positions => null()
        class(Abstract_Particles_Dipolar_Moments), pointer :: dipolar_moments => null()
        class(Abstract_Moved_Positions), pointer :: moved => null()
        class(Abstract_Visitable_Cells), pointer :: intra_cells => null(), inter_cells => null()
        class(Abstract_Pair_Potential), pointer :: wall_pair
        class(Abstract_Ewald_Real_Pair), pointer :: ewald_real_pair => null()
        class(Abstract_Ewald_Real_Particles), pointer :: ewald_real => null()
    end type Move_Candidate

    type, abstract, public :: Abstract_One_Particle_Move
    private
        class(Abstract_Temperature), pointer :: temperature => null()
        class(Abstract_External_Field), pointer :: field => null()
        class(Abstract_Walls_Potential), pointer :: walls => null()
        class(Abstract_Ewald_Real_Pair), pointer :: ewald_inter_real_pair => null()
        type(Move_Candidate) :: candidates(num_components)
        type(Concrete_Change_Counters), pointer :: move_counters(:) => null()
        type(Concrete_Particles_Energy), pointer :: particles_energies(:) => null()
        type(Concrete_Inter_Energy), pointer :: inter_energy => null()
    contains
        procedure :: construct => Abstract_One_Particle_Move_construct
        generic :: set => set_components, set_short_potentials, set_ewalds, set_observables
        procedure :: destroy => Abstract_One_Particle_Move_destroy
        procedure :: try => Abstract_One_Particle_Move_try
        procedure, private :: set_components => Abstract_One_Particle_Move_set_components
        procedure, private :: set_short_potentials => &
            Abstract_One_Particle_Move_set_short_potentials
        procedure, private :: set_ewalds => Abstract_One_Particle_Move_set_ewalds
        procedure, private :: set_observables => Abstract_One_Particle_Move_set_observables
        procedure, private :: test_metropolis => Abstract_One_Particle_Move_test_metropolis
        procedure, private :: nullify_candidate => Abstract_One_Particle_Move_nullify_candidate
        procedure, private, nopass :: select_actor_and_spectator => &
            Abstract_One_Particle_Move_select_actor_and_spectator
        generic, private :: set_actor => set_actor_wall_pair, set_actor_positions_and_moved, &
            set_actor_dipolar_moments, set_actor_cells, set_actor_ewald
        procedure, private :: set_actor_wall_pair => Abstract_One_Particle_Move_set_actor_wall_pair
        procedure, private :: set_actor_positions_and_moved => &
            Abstract_One_Particle_Move_set_actor_positions_and_moved
        procedure, private :: set_actor_dipolar_moments => &
            Abstract_One_Particle_Move_set_actor_dipolar_moments
        procedure, private :: set_actor_cells => Abstract_One_Particle_Move_set_actor_cells
        procedure, private :: set_actor_ewald => Abstract_One_Particle_Move_set_actor_ewald
        procedure, private :: set_spectator => Abstract_One_Particle_Move_set_spectator
    end type Abstract_One_Particle_Move

    type, extends(Abstract_One_Particle_Move), public :: Two_Candidates_One_Particle_Move

    end type Two_Candidates_One_Particle_Move

    type, extends(Abstract_One_Particle_Move), public :: First_Candidate_One_Particle_Move
    contains
        procedure, private, nopass :: select_actor_and_spectator => &
            First_Candidate_One_Particle_Move_select_actor_and_spectator
    end type First_Candidate_One_Particle_Move

    type, extends(Abstract_One_Particle_Move), public :: Second_Candidate_One_Particle_Move
    contains
        procedure, private, nopass :: select_actor_and_spectator => &
            Second_Candidate_One_Particle_Move_select_actor_and_spectator
    end type Second_Candidate_One_Particle_Move

    type, extends(Abstract_One_Particle_Move), public :: Null_One_Particle_Move
    contains
        procedure :: construct => Null_One_Particle_Move_construct
        procedure :: destroy => Null_One_Particle_Move_destroy
        procedure, private :: set_components => Null_One_Particle_Move_set_components
        procedure, private :: set_short_potentials => Null_One_Particle_Move_set_short_potentials
        procedure, private :: set_ewalds => Null_One_Particle_Move_set_ewalds
        procedure, private :: set_observables => Null_One_Particle_Move_set_observables
        procedure :: try => Null_One_Particle_Move_try
    end type Null_One_Particle_Move

contains

!implementation Abstract_One_Particle_Move

    subroutine Abstract_One_Particle_Move_construct(this, environment, moved_1, moved_2)
        class(Abstract_One_Particle_Move), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Moved_Positions), target, intent(in) :: moved_1, moved_2

        this%temperature => environment%temperature
        this%field => environment%external_field
        this%walls => environment%walls_potential
        this%candidates(1)%moved => moved_1
        this%candidates(2)%moved => moved_2
    end subroutine Abstract_One_Particle_Move_construct

    subroutine Abstract_One_Particle_Move_destroy(this)
        class(Abstract_One_Particle_Move), intent(inout) :: this

        this%move_counters => null()
        this%particles_energies => null()
        call this%nullify_candidate(2)
        call this%nullify_candidate(1)
        this%ewald_inter_real_pair => null()
        this%walls => null()
        this%field => null()
        this%temperature => null()
    end subroutine Abstract_One_Particle_Move_destroy

    subroutine Abstract_One_Particle_Move_nullify_candidate(this, i_candidate)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        integer, intent(in) :: i_candidate

        call check_in_range("Abstract_One_Particle_Move_nullify_candidate", num_components, &
            "i_candidate", i_candidate)
        this%candidates(i_candidate)%ewald_real => null()
        this%candidates(i_candidate)%ewald_real_pair => null()
        this%candidates(i_candidate)%inter_cells => null()
        this%candidates(i_candidate)%intra_cells => null()
        this%candidates(i_candidate)%moved => null()
        this%candidates(i_candidate)%dipolar_moments => null()
        this%candidates(i_candidate)%positions => null()
    end subroutine Abstract_One_Particle_Move_nullify_candidate

    subroutine Abstract_One_Particle_Move_set_components(this, components)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(num_components)

        this%candidates(1)%positions => components(1)%positions
        this%candidates(1)%dipolar_moments => components(1)%dipolar_moments
        this%candidates(2)%positions => components(2)%positions
        this%candidates(2)%dipolar_moments => components(2)%dipolar_moments
    end subroutine Abstract_One_Particle_Move_set_components

    subroutine Abstract_One_Particle_Move_set_short_potentials(this, short_potentials)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Mixture_Short_Potentials_Wrapper), target, intent(in) :: short_potentials

        this%candidates(1)%intra_cells => short_potentials%intras(1)%cells
        this%candidates(1)%inter_cells => short_potentials%inters(1)%cells
        this%candidates(1)%wall_pair => short_potentials%intras(1)%wall_pair
        this%candidates(2)%intra_cells => short_potentials%intras(2)%cells
        this%candidates(2)%inter_cells => short_potentials%inters(2)%cells
        this%candidates(2)%wall_pair => short_potentials%intras(2)%wall_pair
    end subroutine Abstract_One_Particle_Move_set_short_potentials

    subroutine Abstract_One_Particle_Move_set_ewalds(this, ewalds)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds

        this%candidates(1)%ewald_real_pair => ewalds%intras(1)%real_pair
        this%candidates(1)%ewald_real => ewalds%intras(1)%real_particles
        this%candidates(2)%ewald_real_pair => ewalds%intras(2)%real_pair
        this%candidates(2)%ewald_real => ewalds%intras(2)%real_particles
        this%ewald_inter_real_pair => ewalds%inter%real_pair
    end subroutine Abstract_One_Particle_Move_set_ewalds

    subroutine Abstract_One_Particle_Move_set_observables(this, move_counters, particles_energies, &
        inter_energy)
        class(Abstract_One_Particle_Move), intent(inout) :: this
        type(Concrete_Change_Counters), target, intent(in) :: move_counters(:)
        type(Concrete_Particles_Energy), target, intent(in) :: particles_energies(:)
        type(Concrete_Inter_Energy), target, intent(in) :: inter_energy

        if (size(move_counters) /= num_components) then
            call error_exit("Abstract_One_Particle_Move: "//&
                "move_counters doesn't have the right size.")
        end if
        this%move_counters => move_counters
        if (size(particles_energies) /= num_components) then
            call error_exit("Abstract_One_Particle_Move: "//&
                "particles_energies doesn't have the right size.")
        end if
        this%particles_energies => particles_energies
        this%inter_energy => inter_energy
    end subroutine Abstract_One_Particle_Move_set_observables

    subroutine Abstract_One_Particle_Move_try(this)
        class(Abstract_One_Particle_Move), intent(in) :: this

        integer :: i_actor, i_spectator
        logical :: success
        type(Concrete_Particles_Energy) :: actor_energy_difference
        type(Concrete_Inter_Energy) :: inter_energy_difference

        call this%select_actor_and_spectator(i_actor, i_spectator)
        this%move_counters(i_actor)%num_hits = this%move_counters(i_actor)%num_hits + 1
        call this%test_metropolis(success, actor_energy_difference, inter_energy_difference, &
            i_actor, i_spectator)
        if (success) then
            this%particles_energies(i_actor) = this%particles_energies(i_actor) + &
                actor_energy_difference
            this%inter_energy = this%inter_energy + inter_energy_difference
            this%move_counters(i_actor)%num_success = this%move_counters(i_actor)%num_success + 1
        end if
    end subroutine Abstract_One_Particle_Move_try

    subroutine Abstract_One_Particle_Move_test_metropolis(this, success, actor_energy_difference, &
        inter_energy_difference, i_actor, i_spectator)
        class(Abstract_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Particles_Energy), intent(out) :: actor_energy_difference
        type(Concrete_Inter_Energy), intent(out) :: inter_energy_difference
        integer, intent(in) :: i_actor, i_spectator

        class(Abstract_Pair_Potential), pointer :: actor_wall_pair
        class(Abstract_Particles_Positions), pointer :: actor_positions
        class(Abstract_Particles_Dipolar_Moments), pointer :: actor_dipolar_moments
        class(Abstract_Moved_Positions), pointer :: actor_moved
        class(Abstract_Visitable_Cells), pointer :: actor_cells, actor_inter_cells, spectator_cells
        class(Abstract_Ewald_Real_Pair), pointer :: actor_ewald_real_pair
        class(Abstract_Ewald_Real_Particles), pointer :: actor_ewald_real, spectator_ewald_real
        integer :: spectator_num_positions
        type(Concrete_Particle) :: old, new
        type(Concrete_Particles_Energy) :: new_energy, old_energy
        type(Concrete_Inter_Energy) :: inter_new_energy, inter_old_energy
        real(DP) :: energy_difference
        logical :: overlap
        real(DP) :: rand

        call this%set_actor(actor_wall_pair, i_actor)
        call this%set_actor(actor_positions, actor_moved, i_actor)
        call this%set_actor(actor_dipolar_moments, i_actor)
        call this%set_actor(actor_cells, actor_inter_cells, i_actor)
        call this%set_actor(actor_ewald_real_pair, actor_ewald_real, i_actor)
        call this%set_spectator(spectator_num_positions, spectator_cells, spectator_ewald_real, &
            i_spectator)

        old%i = random_integer(actor_positions%get_num())
        old%position = actor_positions%get(old%i)
        new%i = old%i
        new%position = actor_moved%get(new%i)

        success = .false.
        call this%walls%visit(overlap, new_energy%walls, new%position, actor_wall_pair)
        if (overlap) return
        if (actor_positions%get_num() > spectator_num_positions) then
            new%same_type = .true.
            call actor_cells%visit(overlap, new_energy%short, new)
            if (overlap) return
            new%same_type = .false.
            call spectator_cells%visit(overlap, inter_new_energy%short, new)
        else
            new%same_type = .false.
            call spectator_cells%visit(overlap, inter_new_energy%short, new)
            if (overlap) return
            new%same_type = .true.
            call actor_cells%visit(overlap, new_energy%short, new)
        end if
        if (overlap) return
        call this%walls%visit(overlap, old_energy%walls, old%position, actor_wall_pair)
        old%same_type = .true.
        call actor_cells%visit(overlap, old_energy%short, old)
        old%same_type = .false.
        call spectator_cells%visit(overlap, inter_old_energy%short, old)

        new%dipolar_moment = actor_dipolar_moments%get(new%i)
        new%same_type = .true.
        call actor_ewald_real%visit(new_energy%long, new, actor_ewald_real_pair)
        new%same_type = .false.
        call spectator_ewald_real%visit(inter_new_energy%long, new, this%ewald_inter_real_pair)
        old%dipolar_moment = actor_dipolar_moments%get(old%i)
        old%same_type = .true.
        call actor_ewald_real%visit(old_energy%long, old, actor_ewald_real_pair)
        old%same_type = .false.
        call spectator_ewald_real%visit(inter_old_energy%long, old, this%ewald_inter_real_pair)

        actor_energy_difference = new_energy - old_energy
        inter_energy_difference = inter_new_energy - inter_old_energy
        energy_difference = particle_energy_sum(actor_energy_difference) + &
            inter_energy_sum(inter_energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference/this%temperature%get())) then
            call actor_positions%set(new%i, new%position)
            call actor_cells%move(old, new)
            call actor_inter_cells%move(old, new)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Move_test_metropolis

    subroutine Abstract_One_Particle_Move_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = random_integer(num_components)
        i_spectator = mod(i_actor, num_components) + 1
    end subroutine Abstract_One_Particle_Move_select_actor_and_spectator

    subroutine Abstract_One_Particle_Move_set_actor_wall_pair(this, actor_wall_pair, i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        class(Abstract_Pair_Potential), pointer, intent(out) :: actor_wall_pair
        integer, intent(in) :: i_actor

        actor_wall_pair => this%candidates(i_actor)%wall_pair
    end subroutine Abstract_One_Particle_Move_set_actor_wall_pair

    subroutine Abstract_One_Particle_Move_set_actor_positions_and_moved(this, actor_positions, &
        actor_moved, i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        class(Abstract_Particles_Positions), pointer, intent(out) :: actor_positions
        class(Abstract_Moved_Positions), pointer, intent(out) :: actor_moved
        integer, intent(in) :: i_actor

        actor_positions => this%candidates(i_actor)%positions
        actor_moved => this%candidates(i_actor)%moved
    end subroutine Abstract_One_Particle_Move_set_actor_positions_and_moved

    subroutine Abstract_One_Particle_Move_set_actor_dipolar_moments(this, actor_dipolar_moments, &
        i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        class(Abstract_Particles_Dipolar_Moments), pointer, intent(out) :: actor_dipolar_moments
        integer, intent(in) :: i_actor

        actor_dipolar_moments => this%candidates(i_actor)%dipolar_moments
    end subroutine Abstract_One_Particle_Move_set_actor_dipolar_moments

    subroutine Abstract_One_Particle_Move_set_actor_cells(this, actor_cells, actor_inter_cells, &
        i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        class(Abstract_Visitable_Cells), pointer, intent(out) :: actor_cells, actor_inter_cells
        integer, intent(in) :: i_actor

        actor_cells => this%candidates(i_actor)%intra_cells
        actor_inter_cells => this%candidates(i_actor)%inter_cells
    end subroutine Abstract_One_Particle_Move_set_actor_cells

    subroutine Abstract_One_Particle_Move_set_actor_ewald(this, actor_ewald_real_pair, &
        actor_ewald_real, i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        class(Abstract_Ewald_Real_Pair), pointer, intent(out) :: actor_ewald_real_pair
        class(Abstract_Ewald_Real_Particles), pointer, intent(out) :: actor_ewald_real
        integer, intent(in) :: i_actor

        actor_ewald_real_pair => this%candidates(i_actor)%ewald_real_pair
        actor_ewald_real => this%candidates(i_actor)%ewald_real
    end subroutine Abstract_One_Particle_Move_set_actor_ewald

    subroutine Abstract_One_Particle_Move_set_spectator(this, spectator_num_positions, &
        spectator_cells, spectator_ewald_real, i_spectator)
        class(Abstract_One_Particle_Move), intent(in) :: this
        integer, intent(out) :: spectator_num_positions
        class(Abstract_Visitable_Cells), pointer, intent(out) :: spectator_cells
        class(Abstract_Ewald_Real_Particles), pointer, intent(out) :: spectator_ewald_real
        integer, intent(in) :: i_spectator

        spectator_num_positions = this%candidates(i_spectator)%positions%get_num()
        spectator_cells => this%candidates(i_spectator)%inter_cells
        spectator_ewald_real => this%candidates(i_spectator)%ewald_real
    end subroutine Abstract_One_Particle_Move_set_spectator

!end implementation Abstract_One_Particle_Move

!implementation First_Candidate_One_Particle_Move

    subroutine First_Candidate_One_Particle_Move_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 1
        i_spectator = 2
    end subroutine First_Candidate_One_Particle_Move_select_actor_and_spectator

!end implementation First_Candidate_One_Particle_Move

!implementation Second_Candidate_One_Particle_Move

    subroutine Second_Candidate_One_Particle_Move_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 2
        i_spectator = 1
    end subroutine Second_Candidate_One_Particle_Move_select_actor_and_spectator

!end implementation Second_Candidate_One_Particle_Move

!implementation Null_One_Particle_Move

    subroutine Null_One_Particle_Move_construct(this, environment, moved_1, moved_2)
        class(Null_One_Particle_Move), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Moved_Positions), target, intent(in) :: moved_1, moved_2
    end subroutine Null_One_Particle_Move_construct

    subroutine Null_One_Particle_Move_destroy(this)
        class(Null_One_Particle_Move), intent(inout) :: this
    end subroutine Null_One_Particle_Move_destroy

    subroutine Null_One_Particle_Move_set_components(this, components)
        class(Null_One_Particle_Move), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(num_components)
    end subroutine Null_One_Particle_Move_set_components

    subroutine Null_One_Particle_Move_set_short_potentials(this, short_potentials)
        class(Null_One_Particle_Move), intent(inout) :: this
        type(Mixture_Short_Potentials_Wrapper), target, intent(in) :: short_potentials
    end subroutine Null_One_Particle_Move_set_short_potentials

    subroutine Null_One_Particle_Move_set_ewalds(this, ewalds)
        class(Null_One_Particle_Move), intent(inout) :: this
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds
    end subroutine Null_One_Particle_Move_set_ewalds

    subroutine Null_One_Particle_Move_set_observables(this, move_counters, particles_energies, &
        inter_energy)
        class(Null_One_Particle_Move), intent(inout) :: this
        type(Concrete_Change_Counters), target, intent(in) :: move_counters(:)
        type(Concrete_Particles_Energy), target, intent(in) :: particles_energies(:)
        type(Concrete_Inter_Energy), target, intent(in) :: inter_energy
    end subroutine Null_One_Particle_Move_set_observables

    subroutine Null_One_Particle_Move_try(this)
        class(Null_One_Particle_Move), intent(in) :: this
    end subroutine Null_One_Particle_Move_try

!end implementation Null_One_Particle_Move

end module class_one_particle_move
