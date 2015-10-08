module class_one_particle_rotation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components
use procedures_errors, only: error_exit
use procedures_checks, only: check_in_range
use procedures_random, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use class_temperature, only: Abstract_Temperature
use class_external_field, only: Abstract_External_Field
use class_particles_moment_norm, only: Abstract_Particles_Moment_Norm
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_orientations, only: Abstract_Particles_Orientations
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use types_particles_wrapper, only: Particles_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_rotated_orientations, only: Abstract_Rotated_Orientations
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

    type :: Rotation_Candidate
        class(Abstract_Particles_Moment_Norm), pointer :: moment_norm => null()
        class(Abstract_Particles_Positions), pointer :: positions => null()
        class(Abstract_Particles_Orientations), pointer :: orientations => null()
        class(Abstract_Particles_Dipolar_Moments), pointer :: dipolar_moments => null()
        class(Abstract_Rotated_Orientations), pointer :: rotated => null()
        class(Abstract_Ewald_Real_Pair), pointer :: ewald_real_pair => null()
        class(Abstract_Ewald_Real_Particles), pointer :: ewald_real => null()
    end type Rotation_Candidate

    type, abstract, public :: Abstract_One_Particle_Rotation
    private
        class(Abstract_Temperature), pointer :: temperature => null()
        class(Abstract_External_Field), pointer :: field => null()
        class(Abstract_Ewald_Real_Pair), pointer :: ewald_inter_real_pair => null()
        type(Rotation_Candidate) :: candidates(num_components)
        type(Concrete_Change_Counters), pointer :: rotation_counters(:) => null()
        type(Concrete_Particles_Energy), pointer :: particles_energies(:) => null()
        type(Concrete_Inter_Energy), pointer :: inter_energy => null()
    contains
        procedure :: construct => Abstract_One_Particle_Rotation_construct
        generic :: set => set_components, set_ewalds, set_observables
        procedure :: destroy => Abstract_One_Particle_Rotation_destroy
        procedure :: try => Abstract_One_Particle_Rotation_try
        procedure, private :: set_components => Abstract_One_Particle_Rotation_set_components
        procedure, private :: set_ewalds => Abstract_One_Particle_Rotation_set_ewalds
        procedure, private :: set_observables => Abstract_One_Particle_Rotation_set_observables
        procedure, private :: test_metropolis => Abstract_One_Particle_Rotation_test_metropolis
        procedure, private :: nullify_candidate => Abstract_One_Particle_Rotation_nullify_candidate
        procedure, private, nopass :: select_actor_and_spectator => &
            Abstract_One_Particle_Rotation_select_actor_and_spectator
        procedure, private :: visit_long => Abstract_One_Particle_Rotation_visit_long
    end type Abstract_One_Particle_Rotation

    type, extends(Abstract_One_Particle_Rotation), public :: Two_Candidates_Rotation

    end type Two_Candidates_Rotation

    type, extends(Abstract_One_Particle_Rotation), public :: First_Candidate_Rotation
    contains
        procedure, private, nopass :: select_actor_and_spectator => &
            First_Rotation_select_actor_and_spectator
    end type First_Candidate_Rotation

    type, extends(Abstract_One_Particle_Rotation), public :: Second_Candidate_Rotation
    contains
        procedure, private, nopass :: select_actor_and_spectator => &
            Second_Candidate_Rotation_select_actor_and_spectator
    end type Second_Candidate_Rotation

    type, extends(Abstract_One_Particle_Rotation), public :: Null_One_Particle_Rotation
    contains
        procedure :: construct => Null_One_Particle_Rotation_construct
        procedure :: destroy => Null_One_Particle_Rotation_destroy
        procedure :: try => Null_One_Particle_Rotation_try
        procedure, private :: set_components => Null_One_Particle_Rotation_set_components
        procedure, private :: set_ewalds => Null_One_Particle_Rotation_set_ewalds
        procedure, private :: set_observables => Null_One_Particle_Rotation_set_observables
    end type Null_One_Particle_Rotation

contains

!implementation Abstract_One_Particle_Rotation

    subroutine Abstract_One_Particle_Rotation_construct(this, environment, rotated_1, rotated_2)
        class(Abstract_One_Particle_Rotation), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Rotated_Orientations), target, intent(in) :: rotated_1, rotated_2

        this%temperature => environment%temperature
        this%field => environment%external_field
        this%candidates(1)%rotated => rotated_1
        this%candidates(2)%rotated => rotated_2
    end subroutine Abstract_One_Particle_Rotation_construct

    subroutine Abstract_One_Particle_Rotation_destroy(this)
        class(Abstract_One_Particle_Rotation), intent(inout) :: this

        this%rotation_counters => null()
        this%particles_energies => null()
        call this%nullify_candidate(2)
        call this%nullify_candidate(1)
        this%ewald_inter_real_pair => null()
        this%field => null()
        this%temperature => null()
    end subroutine Abstract_One_Particle_Rotation_destroy

    subroutine Abstract_One_Particle_Rotation_nullify_candidate(this, i_candidate)
        class(Abstract_One_Particle_Rotation), intent(inout) :: this
        integer, intent(in) :: i_candidate

        call check_in_range("Abstract_One_Particle_Rotation_nullify_candidate", num_components, &
            "i_candidate", i_candidate)
        this%candidates(i_candidate)%ewald_real => null()
        this%candidates(i_candidate)%ewald_real_pair => null()
        this%candidates(i_candidate)%rotated => null()
        this%candidates(i_candidate)%dipolar_moments => null()
        this%candidates(i_candidate)%orientations => null()
        this%candidates(i_candidate)%positions => null()
        this%candidates(i_candidate)%moment_norm => null()
    end subroutine Abstract_One_Particle_Rotation_nullify_candidate

    subroutine Abstract_One_Particle_Rotation_set_components(this, components)
        class(Abstract_One_Particle_Rotation), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(num_components)

        integer :: i_candidate
        do i_candidate = 1, num_components
            this%candidates(i_candidate)%moment_norm => components(i_candidate)%moment_norm
            this%candidates(i_candidate)%positions => components(i_candidate)%positions
            this%candidates(i_candidate)%orientations => components(i_candidate)%orientations
            this%candidates(i_candidate)%dipolar_moments => components(i_candidate)%dipolar_moments
        end do
    end subroutine Abstract_One_Particle_Rotation_set_components

    subroutine Abstract_One_Particle_Rotation_set_ewalds(this, ewalds)
        class(Abstract_One_Particle_Rotation), intent(inout) :: this
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds

        integer :: i_candidate
        do i_candidate = 1, num_components
            this%candidates(i_candidate)%ewald_real_pair => ewalds%intras(i_candidate)%real_pair
            this%candidates(i_candidate)%ewald_real => ewalds%intras(i_candidate)%real_particles
        end do
        this%ewald_inter_real_pair => ewalds%inter%real_pair
    end subroutine Abstract_One_Particle_Rotation_set_ewalds

    subroutine Abstract_One_Particle_Rotation_set_observables(this, rotation_counters, &
        particles_energies, inter_energy)
        class(Abstract_One_Particle_Rotation), intent(inout) :: this
        type(Concrete_Change_Counters), target, intent(in) :: rotation_counters(:)
        type(Concrete_Particles_Energy), target, intent(in) :: particles_energies(:)
        type(Concrete_Inter_Energy), target, intent(in) :: inter_energy

        if (size(rotation_counters) /= num_components) then
            call error_exit("Abstract_One_Particle_Rotation: "//&
                "rotation_counters doesn't have the right size.")
        end if
        this%rotation_counters => rotation_counters
        if (size(particles_energies) /= num_components) then
            call error_exit("Abstract_One_Particle_Rotation: "//&
                "particles_energies doesn't have the right size.")
        end if
        this%particles_energies => particles_energies
        this%inter_energy => inter_energy
    end subroutine Abstract_One_Particle_Rotation_set_observables

    subroutine Abstract_One_Particle_Rotation_try(this)
        class(Abstract_One_Particle_Rotation), intent(in) :: this

        integer :: i_actor, i_spectator
        logical :: success
        type(Concrete_Particles_Energy) :: actor_energy_difference
        type(Concrete_Inter_Energy) :: inter_energy_difference

        call this%select_actor_and_spectator(i_actor, i_spectator)
        this%rotation_counters(i_actor)%num_hits = this%rotation_counters(i_actor)%num_hits + 1
        call this%test_metropolis(success, actor_energy_difference, inter_energy_difference, &
            i_actor, i_spectator)
        if (success) then
            this%particles_energies(i_actor) = this%particles_energies(i_actor) + &
                actor_energy_difference
            this%inter_energy = this%inter_energy + inter_energy_difference
            this%rotation_counters(i_actor)%num_success = &
                this%rotation_counters(i_actor)%num_success + 1
        end if
    end subroutine Abstract_One_Particle_Rotation_try

    subroutine Abstract_One_Particle_Rotation_test_metropolis(this, success, &
        actor_energy_difference, inter_energy_difference, i_actor, i_spectator)
        class(Abstract_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Particles_Energy), intent(out) :: actor_energy_difference
        type(Concrete_Inter_Energy), intent(out) :: inter_energy_difference
        integer, intent(in) :: i_actor, i_spectator

        type(Concrete_Temporary_Particle) :: new, old
        type(Concrete_Mixture_Energy) :: new_energy, old_energy
        real(DP) :: energy_difference
        real(DP) :: rand

        old%i = random_integer(this%candidates(i_actor)%orientations%get_num())
        old%position = this%candidates(i_actor)%positions%get(old%i)
        old%dipolar_moment = this%candidates(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%candidates(i_actor)%rotated%get(new%i)
        new%dipolar_moment = this%candidates(i_actor)%moment_norm%get() * new%orientation

        success = .false.
        call this%visit_long(new_energy, old_energy, new, old, i_actor, i_spectator)

        actor_energy_difference = new_energy%intras(i_actor) - old_energy%intras(i_actor)
        inter_energy_difference = new_energy%inter - old_energy%inter
        energy_difference = particle_energy_sum(actor_energy_difference) + &
            inter_energy_sum(inter_energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference/this%temperature%get())) then
            call this%candidates(i_actor)%orientations%set(new%i, new%orientation)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Rotation_test_metropolis

    subroutine Abstract_One_Particle_Rotation_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = random_integer(num_components)
        i_spectator = mod(i_actor, num_components) + 1
    end subroutine Abstract_One_Particle_Rotation_select_actor_and_spectator

    subroutine Abstract_One_Particle_Rotation_visit_long(this, new_energy, old_energy, new, old, &
        i_actor, i_spectator)
        class(Abstract_One_Particle_Rotation), intent(in) :: this
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
    end subroutine Abstract_One_Particle_Rotation_visit_long

!end implementation Abstract_One_Particle_Rotation

!implementation First_Candidate_Rotation

    subroutine First_Rotation_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 1
        i_spectator = 2
    end subroutine First_Rotation_select_actor_and_spectator

!end implementation First_Candidate_Rotation

!implementation Second_Candidate_Rotation

    subroutine Second_Candidate_Rotation_select_actor_and_spectator(i_actor, i_spectator)
        integer, intent(out) :: i_actor, i_spectator

        i_actor = 2
        i_spectator = 1
    end subroutine Second_Candidate_Rotation_select_actor_and_spectator

!end implementation Second_Candidate_Rotation

!implementation Null_One_Particle_Rotation

    subroutine Null_One_Particle_Rotation_construct(this, environment, rotated_1, rotated_2)
        class(Null_One_Particle_Rotation), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Rotated_Orientations), target, intent(in) :: rotated_1, rotated_2
    end subroutine Null_One_Particle_Rotation_construct

    subroutine Null_One_Particle_Rotation_destroy(this)
        class(Null_One_Particle_Rotation), intent(inout) :: this
    end subroutine Null_One_Particle_Rotation_destroy

    subroutine Null_One_Particle_Rotation_set_components(this, components)
        class(Null_One_Particle_Rotation), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(num_components)
    end subroutine Null_One_Particle_Rotation_set_components

    subroutine Null_One_Particle_Rotation_set_ewalds(this, ewalds)
        class(Null_One_Particle_Rotation), intent(inout) :: this
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds
    end subroutine Null_One_Particle_Rotation_set_ewalds

    subroutine Null_One_Particle_Rotation_set_observables(this, rotation_counters, &
        particles_energies, inter_energy)
        class(Null_One_Particle_Rotation), intent(inout) :: this
        type(Concrete_Change_Counters), target, intent(in) :: rotation_counters(:)
        type(Concrete_Particles_Energy), target, intent(in) :: particles_energies(:)
        type(Concrete_Inter_Energy), target, intent(in) :: inter_energy
    end subroutine Null_One_Particle_Rotation_set_observables

    subroutine Null_One_Particle_Rotation_try(this)
        class(Null_One_Particle_Rotation), intent(in) :: this
    end subroutine Null_One_Particle_Rotation_try

!end implementation Null_One_Particle_Rotation

end module class_one_particle_rotation
