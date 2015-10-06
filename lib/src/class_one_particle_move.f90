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
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use types_particles_wrapper, only: Particles_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_moved_positions, only: Abstract_Moved_Positions
use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_cells, only: Abstract_Visitable_Cells
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_ewald_real_particles, only: Abstract_Ewald_Real_Particles
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use module_particles_energy, only: Concrete_Particles_Energy, Concrete_Inter_Energy, &
    Concrete_Long_Energy, Concrete_Mixture_Energy, &
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
        procedure, private :: visit_short => Abstract_One_Particle_Move_visit_short
        procedure, private :: visit_long => Abstract_One_Particle_Move_visit_long
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

        type(Concrete_Temporary_Particle) :: new, old
        type(Concrete_Mixture_Energy) :: new_energy, old_energy
        real(DP) :: energy_difference
        logical :: overlap
        real(DP) :: rand

        old%i = random_integer(this%candidates(i_actor)%positions%get_num())
        old%position = this%candidates(i_actor)%positions%get(old%i)
        old%dipolar_moment = this%candidates(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = this%candidates(i_actor)%moved%get(new%i)
        new%dipolar_moment = old%dipolar_moment

        success = .false.
        call this%visit_short(overlap, new_energy, old_energy, new, old, i_actor, i_spectator)
        if (overlap) return
        call this%visit_long(new_energy, old_energy, new, old, i_actor, i_spectator)

        actor_energy_difference = new_energy%intras(i_actor) - old_energy%intras(i_actor)
        inter_energy_difference = new_energy%inter - old_energy%inter
        energy_difference = particle_energy_sum(actor_energy_difference) + &
            inter_energy_sum(inter_energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference/this%temperature%get())) then
            call this%candidates(i_actor)%positions%set(new%i, new%position)
            call this%candidates(i_actor)%intra_cells%move(old, new)
            call this%candidates(i_actor)%inter_cells%move(old, new)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Move_test_metropolis

    subroutine Abstract_One_Particle_Move_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = random_integer(num_components)
        i_spectator = mod(i_actor, num_components) + 1
    end subroutine Abstract_One_Particle_Move_select_actor_and_spectator

    subroutine Abstract_One_Particle_Move_visit_short(this, overlap, new_energy, old_energy,  new, &
        old, i_actor, i_spectator)
        class(Abstract_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Mixture_Energy), intent(out) :: new_energy, old_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor, i_spectator

        associate(actor_new_energy => new_energy%intras(i_actor), &
            actor_old_energy => old_energy%intras(i_actor), &
            actor_wall_pair => this%candidates(i_actor)%wall_pair, &
            actor_num_positions => this%candidates(i_actor)%positions%get_num(), &
            spectator_num_positions => this%candidates(i_spectator)%positions%get_num(), &
            actor_cells => this%candidates(i_actor)%intra_cells, &
            spectator_cells => this%candidates(i_spectator)%inter_cells)

            call this%walls%visit(overlap, actor_new_energy%walls, new%position, actor_wall_pair)
            if (overlap) return
            if (actor_num_positions > spectator_num_positions) then
                call actor_cells%visit(overlap, actor_new_energy%short, new, same_type=.true.)
                if (overlap) return
                call spectator_cells%visit(overlap, new_energy%inter%short, new, same_type=.false.)
            else
                call spectator_cells%visit(overlap, new_energy%inter%short, new, same_type=.false.)
                if (overlap) return
                call actor_cells%visit(overlap, actor_new_energy%short, new, same_type=.true.)
            end if
            if (overlap) return
            call this%walls%visit(overlap, actor_old_energy%walls, old%position, actor_wall_pair)
            call actor_cells%visit(overlap, actor_old_energy%short, old, same_type=.true.)
            call spectator_cells%visit(overlap, old_energy%inter%short, old, same_type=.false.)
        end associate
    end subroutine Abstract_One_Particle_Move_visit_short

    subroutine Abstract_One_Particle_Move_visit_long(this, new_energy, old_energy, new, old, &
        i_actor, i_spectator)
        class(Abstract_One_Particle_Move), intent(in) :: this
        type(Concrete_Mixture_Energy), intent(inout) :: new_energy, old_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor, i_spectator

        type(Concrete_Long_Energy) :: new_long, old_long, inter_new_long, inter_old_long

        associate(actor_ewald_real_pair => this%candidates(i_actor)%ewald_real_pair, &
            actor_ewald_real => this%candidates(i_actor)%ewald_real, &
            spectator_ewald_real => this%candidates(i_spectator)%ewald_real)

            call actor_ewald_real%visit(new_long%real, new, actor_ewald_real_pair, &
                same_type=.true.)
            call spectator_ewald_real%visit(inter_new_long%real, new, &
                this%ewald_inter_real_pair, same_type=.false.)
            call actor_ewald_real%visit(old_long%real, old, actor_ewald_real_pair, &
                same_type=.true.)
            call spectator_ewald_real%visit(inter_old_long%real, old, &
                this%ewald_inter_real_pair, same_type=.false.)
        end associate
        new_energy%intras(i_actor)%long = new_long%real
        old_energy%intras(i_actor)%long = old_long%real
        new_energy%inter%long = inter_new_long%real
        old_energy%inter%long = inter_old_long%real
    end subroutine Abstract_One_Particle_Move_visit_long

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
