module class_one_particle_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_components
use procedures_errors, only: error_exit
use procedures_checks, only: check_in_range
use procedures_random, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_particles_wrapper, only: Particles_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_changed_coordinates, only: Abstract_Changed_Coordinates
use types_short_potential_wrapper, only: Mixture_Short_Potentials_Wrapper
use types_ewald_wrapper, only: Mixture_Ewald_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler
use module_changes_success, only: Concrete_Change_Counter
use module_particles_energy, only: Concrete_Particles_Energy, Concrete_Inter_Energy, &
    Concrete_Long_Energy, Concrete_Mixture_Energy, &
    particle_energy_sum => Concrete_Particles_Energy_sum, &
    inter_energy_sum => Concrete_Inter_Energy_sum, operator(+), operator(-)
use types_observables_wrapper, only: Mixture_Observables_Wrapper
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_One_Particle_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Particles_Wrapper), pointer :: components(:) => null()
        class(Abstract_Changed_Coordinates), pointer :: changed_coordinates(:) => null()
        type(Mixture_Short_Potentials_Wrapper), pointer :: short_potentials => null()
        type(Mixture_Ewald_Wrapper), pointer :: ewalds => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
        type(Mixture_Observables_Wrapper), pointer :: observables
    contains
        procedure :: construct => Abstract_One_Particle_Change_construct
        generic :: set => set_candidates, set_observables
        procedure :: destroy => Abstract_One_Particle_Change_destroy
        procedure :: try => Abstract_One_Particle_Change_try
        procedure(Abstract_One_Particle_Change_set_selector), private, deferred :: set_selector
        procedure, private :: set_candidates => Abstract_One_Particle_Change_set_candidates
        procedure, private :: set_observables => Abstract_One_Particle_Change_set_observables
        procedure, private :: test_metropolis => Abstract_One_Particle_Change_test_metropolis
        procedure(Abstract_One_Particle_Change_define_change), private, deferred :: define_change
        procedure, private :: visit_short => Abstract_One_Particle_Change_visit_short
        procedure, private :: visit_long => Abstract_One_Particle_Change_visit_long
        procedure(Abstract_One_Particle_Change_update_actor), private, deferred :: update_actor
        procedure(Abstract_One_Particle_Change_increment_hits), private, deferred :: increment_hits
        procedure(Abstract_One_Particle_Change_increment_success), private, deferred :: &
            increment_success
    end type Abstract_One_Particle_Change

    abstract interface

        subroutine Abstract_One_Particle_Change_set_selector(this)
        import :: Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(inout) :: this
        end subroutine Abstract_One_Particle_Change_set_selector

        subroutine Abstract_One_Particle_Change_define_change(this, new, old, i_actor)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            type(Concrete_Temporary_Particle), intent(out) :: new, old
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_define_change

        subroutine Abstract_One_Particle_Change_update_actor(this, i_actor, new, old)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(out) :: new, old
        end subroutine Abstract_One_Particle_Change_update_actor

        subroutine Abstract_One_Particle_Change_increment_hits(this, i_actor)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(inout) :: this
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_increment_hits

        subroutine Abstract_One_Particle_Change_increment_success(this, i_actor)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(inout) :: this
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_increment_success

    end interface

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Move
    contains
        procedure :: get_num_choices => Concrete_One_Particle_Move_get_num_choices
        procedure, private :: set_selector => Concrete_One_Particle_Move_set_selector
        procedure, private :: define_change => Concrete_One_Particle_Move_define_change
        procedure, private :: update_actor => Concrete_One_Particle_Move_update_actor
        procedure, private :: increment_hits => Concrete_One_Particle_Move_increment_hits
        procedure, private :: increment_success => Concrete_One_Particle_Move_increment_success
    end type Concrete_One_Particle_Move

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Rotation
    contains
        procedure :: get_num_choices => Concrete_One_Particle_Rotation_get_num_choices
        procedure, private :: set_selector => Concrete_One_Particle_Rotation_set_selector
        procedure, private :: visit_short => Concrete_One_Particle_Rotation_visit_short
        procedure, private :: define_change => Concrete_One_Particle_Rotation_define_change
        procedure, private :: update_actor => Concrete_One_Particle_Rotation_update_actor
        procedure, private :: increment_hits => Concrete_One_Particle_Rotation_increment_hits
        procedure, private :: increment_success => Concrete_One_Particle_Rotation_increment_success
    end type Concrete_One_Particle_Rotation

    type, extends(Abstract_One_Particle_Change), public :: Null_One_Particle_Change
    contains
        procedure :: construct => Null_One_Particle_Change_construct
        procedure :: destroy => Null_One_Particle_Change_destroy
        procedure, private :: set_candidates => Null_One_Particle_Change_set_candidates
        procedure, private :: set_observables => Null_One_Particle_Change_set_observables
        procedure :: get_num_choices => Null_One_Particle_Change_get_num_choices
        procedure :: try => Null_One_Particle_Change_try
        procedure, private :: set_selector => Null_One_Particle_Change_set_selector
        procedure, private :: define_change => Null_One_Particle_Change_define_change
        procedure, private :: update_actor => Null_One_Particle_Change_update_actor
        procedure, private :: increment_hits => Null_One_Particle_Change_increment_hits
        procedure, private :: increment_success => Null_One_Particle_Change_increment_success
    end type Null_One_Particle_Change

contains

!implementation Abstract_One_Particle_Change

    subroutine Abstract_One_Particle_Change_construct(this, environment, changed_coordinates, &
        selector)
        class(Abstract_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Changed_Coordinates), target, intent(in) :: &
            changed_coordinates(num_components)
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%environment => environment
        this%changed_coordinates => changed_coordinates
        allocate(this%selector, source=selector)
    end subroutine Abstract_One_Particle_Change_construct

    subroutine Abstract_One_Particle_Change_destroy(this)
        class(Abstract_One_Particle_Change), intent(inout) :: this

        this%observables => null()
        this%ewalds => null()
        this%short_potentials => null()
        this%components => null()
        call this%selector%destroy()
        if (allocated(this%selector)) deallocate(this%selector)
        this%changed_coordinates => null()
        this%environment => null()
    end subroutine Abstract_One_Particle_Change_destroy

    subroutine Abstract_One_Particle_Change_set_candidates(this, components, short_potentials, &
        ewalds)
        class(Abstract_One_Particle_Change), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(:)
        type(Mixture_Short_Potentials_Wrapper), target, intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds

        if (size(components) /= num_components) then
            call error_exit("Abstract_One_Particle_Change: "//&
                "components doesn't have the right size.")
        end if
        this%components => components
        call this%set_selector()
        this%short_potentials => short_potentials
        this%ewalds => ewalds
    end subroutine Abstract_One_Particle_Change_set_candidates

    subroutine Abstract_One_Particle_Change_set_observables(this, observables)
        class(Abstract_One_Particle_Change), intent(inout) :: this
        type(Mixture_Observables_Wrapper), target, intent(in) :: observables

        this%observables => observables
    end subroutine Abstract_One_Particle_Change_set_observables

    subroutine Abstract_One_Particle_Change_try(this)
        class(Abstract_One_Particle_Change), intent(inout) :: this

        integer :: i_actor, i_spectator
        logical :: success
        type(Concrete_Particles_Energy) :: actor_energy_difference
        type(Concrete_Inter_Energy) :: inter_energy_difference

        i_actor = this%selector%get()
        i_spectator = mod(i_actor, num_components) + 1
        call this%increment_hits(i_actor)
        call this%test_metropolis(success, actor_energy_difference, inter_energy_difference, &
            i_actor, i_spectator)
        if (success) then
            this%observables%intras(i_actor)%particles_energy = &
                this%observables%intras(i_actor)%particles_energy + actor_energy_difference
            this%observables%inter_energy = this%observables%inter_energy + inter_energy_difference
            call this%increment_success(i_actor)
        end if
    end subroutine Abstract_One_Particle_Change_try

    subroutine Abstract_One_Particle_Change_test_metropolis(this, success, &
        actor_energy_difference, inter_energy_difference, i_actor, i_spectator)
        class(Abstract_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Particles_Energy), intent(out) :: actor_energy_difference
        type(Concrete_Inter_Energy), intent(out) :: inter_energy_difference
        integer, intent(in) :: i_actor, i_spectator

        type(Concrete_Temporary_Particle) :: new, old
        type(Concrete_Mixture_Energy) :: new_energy, old_energy
        real(DP) :: energy_difference
        logical :: overlap
        real(DP) :: rand

        call this%define_change(new, old, i_actor)

        success = .false.
        call this%visit_short(overlap, new_energy, old_energy, new, old, i_actor, i_spectator)
        if (overlap) return
        call this%visit_long(new_energy, old_energy, new, old, i_actor, i_spectator)

        actor_energy_difference = new_energy%intras(i_actor) - old_energy%intras(i_actor)
        inter_energy_difference = new_energy%inter - old_energy%inter
        energy_difference = particle_energy_sum(actor_energy_difference) + &
            inter_energy_sum(inter_energy_difference)
        call random_number(rand)
        if (rand < exp(-energy_difference/this%environment%temperature%get())) then
            call this%update_actor(i_actor, new, old)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Change_test_metropolis

    subroutine Abstract_One_Particle_Change_visit_short(this, overlap, new_energy, old_energy, &
        new, old, i_actor, i_spectator)
        class(Abstract_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Mixture_Energy), intent(out) :: new_energy, old_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor, i_spectator

        associate(actor_new_energy => new_energy%intras(i_actor), &
            actor_old_energy => old_energy%intras(i_actor), &
            wall_potential => this%environment%walls_potential, &
            actor_wall_pair => this%short_potentials%intras(i_actor)%wall_pair, &
            actor_num_particles => this%components(i_actor)%number%get(), &
            spectator_num_particles => this%components(i_spectator)%number%get(), &
            actor_cells => this%short_potentials%intras(i_actor)%cells, &
            spectator_cells => this%short_potentials%inters(i_spectator)%cells)

            call wall_potential%visit(overlap, actor_new_energy%walls, new%position, &
                actor_wall_pair)
            if (overlap) return
            if (actor_num_particles > spectator_num_particles) then
                call actor_cells%visit(overlap, actor_new_energy%short, new, same_type=.true.)
                if (overlap) return
                call spectator_cells%visit(overlap, new_energy%inter%short, new, same_type=.false.)
            else
                call spectator_cells%visit(overlap, new_energy%inter%short, new, same_type=.false.)
                if (overlap) return
                call actor_cells%visit(overlap, actor_new_energy%short, new, same_type=.true.)
            end if
            if (overlap) return
            call wall_potential%visit(overlap, actor_old_energy%walls, old%position, &
                actor_wall_pair)
            call actor_cells%visit(overlap, actor_old_energy%short, old, same_type=.true.)
            call spectator_cells%visit(overlap, old_energy%inter%short, old, same_type=.false.)
        end associate
    end subroutine Abstract_One_Particle_Change_visit_short

    subroutine Abstract_One_Particle_Change_visit_long(this, new_energy, old_energy, new, old, &
        i_actor, i_spectator)
        class(Abstract_One_Particle_Change), intent(in) :: this
        type(Concrete_Mixture_Energy), intent(inout) :: new_energy, old_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor, i_spectator

        type(Concrete_Long_Energy) :: new_long, old_long, inter_new_long, inter_old_long

        associate(actor_ewald_real_pair => this%ewalds%intras(i_actor)%real_pair, &
            actor_ewald_real => this%ewalds%intras(i_actor)%real_particles, &
            spectator_ewald_real => this%ewalds%intras(i_spectator)%real_particles)

            call actor_ewald_real%visit(new_long%real, new, actor_ewald_real_pair, &
                same_type=.true.)
            call spectator_ewald_real%visit(inter_new_long%real, new, &
                this%ewalds%inter%real_pair, same_type=.false.)
            call actor_ewald_real%visit(old_long%real, old, actor_ewald_real_pair, &
                same_type=.true.)
            call spectator_ewald_real%visit(inter_old_long%real, old, &
                this%ewalds%inter%real_pair, same_type=.false.)
        end associate
        new_energy%intras(i_actor)%long = new_long%real
        old_energy%intras(i_actor)%long = old_long%real
        new_energy%inter%long = inter_new_long%real
        old_energy%inter%long = inter_old_long%real
        ! add effect of field
    end subroutine Abstract_One_Particle_Change_visit_long

!end implementation Abstract_One_Particle_Change

!implementation Concrete_One_Particle_Move

    subroutine Concrete_One_Particle_Move_set_selector(this)
        class(Concrete_One_Particle_Move), intent(inout) :: this

        call this%selector%construct([this%components(1)%positions%get_num(), &
            this%components(2)%positions%get_num()])
    end subroutine Concrete_One_Particle_Move_set_selector

    pure integer function Concrete_One_Particle_Move_get_num_choices(this) result(num_choices)
        class(Concrete_One_Particle_Move), intent(in) :: this

        num_choices = this%components(1)%positions%get_num() + &
            this%components(2)%positions%get_num()
    end function Concrete_One_Particle_Move_get_num_choices

    subroutine Concrete_One_Particle_Move_define_change(this, new, old, i_actor)
        class(Concrete_One_Particle_Move), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        old%i = random_integer(this%components(i_actor)%positions%get_num())
        old%position = this%components(i_actor)%positions%get(old%i)
        old%dipolar_moment = this%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = this%changed_coordinates(i_actor)%get(new%i)
        new%dipolar_moment = old%dipolar_moment
    end subroutine Concrete_One_Particle_Move_define_change

    subroutine Concrete_One_Particle_Move_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(out) :: new, old

        call this%components(i_actor)%positions%set(new%i, new%position)
        call this%short_potentials%intras(i_actor)%cells%move(old, new)
        call this%short_potentials%inters(i_actor)%cells%move(old, new)
    end subroutine Concrete_One_Particle_Move_update_actor

    subroutine Concrete_One_Particle_Move_increment_hits(this, i_actor)
        class(Concrete_One_Particle_Move), intent(inout) :: this
        integer, intent(in) :: i_actor

        associate (actor_move => this%observables%intras(i_actor)%changes_counter%move)
            actor_move%num_hits = actor_move%num_hits + 1
        end associate
    end subroutine Concrete_One_Particle_Move_increment_hits

    subroutine Concrete_One_Particle_Move_increment_success(this, i_actor)
        class(Concrete_One_Particle_Move), intent(inout) :: this
        integer, intent(in) :: i_actor

        associate (actor_move => this%observables%intras(i_actor)%changes_counter%move)
            actor_move%num_success = actor_move%num_success + 1
        end associate
    end subroutine Concrete_One_Particle_Move_increment_success

!end implementation Concrete_One_Particle_Move

!implementation Concrete_One_Particle_Rotation

    subroutine Concrete_One_Particle_Rotation_set_selector(this)
        class(Concrete_One_Particle_Rotation), intent(inout) :: this

        call this%selector%construct([this%components(1)%orientations%get_num(), &
            this%components(2)%orientations%get_num()])
    end subroutine Concrete_One_Particle_Rotation_set_selector

    pure integer function Concrete_One_Particle_Rotation_get_num_choices(this) result(num_choices)
        class(Concrete_One_Particle_Rotation), intent(in) :: this

        num_choices = this%components(1)%orientations%get_num() + &
            this%components(2)%orientations%get_num()
    end function Concrete_One_Particle_Rotation_get_num_choices

    subroutine Concrete_One_Particle_Rotation_visit_short(this, overlap, new_energy, old_energy, &
        new, old, i_actor, i_spectator)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Mixture_Energy), intent(out) :: new_energy, old_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor, i_spectator
        overlap = .false.
        new_energy%intras(:)%short = 0._DP
        new_energy%inter%short = 0._DP
        old_energy%intras(:)%short = 0._DP
        old_energy%inter%short = 0._DP
    end subroutine Concrete_One_Particle_Rotation_visit_short

    subroutine Concrete_One_Particle_Rotation_define_change(this, new, old, i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        old%i = random_integer(this%components(i_actor)%orientations%get_num())
        old%position = this%components(i_actor)%positions%get(old%i)
        old%dipolar_moment = this%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%changed_coordinates(i_actor)%get(new%i)
        new%dipolar_moment = this%components(i_actor)%moment_norm%get() * new%orientation
    end subroutine Concrete_One_Particle_Rotation_define_change

    subroutine Concrete_One_Particle_Rotation_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(out) :: new, old

        call this%components(i_actor)%orientations%set(new%i, new%orientation)
    end subroutine Concrete_One_Particle_Rotation_update_actor

    subroutine Concrete_One_Particle_Rotation_increment_hits(this, i_actor)
        class(Concrete_One_Particle_Rotation), intent(inout) :: this
        integer, intent(in) :: i_actor

        associate (actor_rotation => this%observables%intras(i_actor)%changes_counter%rotation)
            actor_rotation%num_hits = actor_rotation%num_hits + 1
        end associate
    end subroutine Concrete_One_Particle_Rotation_increment_hits

    subroutine Concrete_One_Particle_Rotation_increment_success(this, i_actor)
        class(Concrete_One_Particle_Rotation), intent(inout) :: this
        integer, intent(in) :: i_actor

        associate (actor_rotation => this%observables%intras(i_actor)%changes_counter%rotation)
            actor_rotation%num_success = actor_rotation%num_success + 1
        end associate
    end subroutine Concrete_One_Particle_Rotation_increment_success

!end implementation Concrete_One_Particle_Rotation

!implementation Null_One_Particle_Change

    subroutine Null_One_Particle_Change_construct(this, environment, changed_coordinates, &
        selector)
        class(Null_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Changed_Coordinates), target, intent(in) :: &
            changed_coordinates(num_components)
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_One_Particle_Change_construct

    subroutine Null_One_Particle_Change_destroy(this)
        class(Null_One_Particle_Change), intent(inout) :: this
    end subroutine Null_One_Particle_Change_destroy

    subroutine Null_One_Particle_Change_set_candidates(this, components, short_potentials, &
        ewalds)
        class(Null_One_Particle_Change), intent(inout) :: this
        type(Particles_Wrapper), target, intent(in) :: components(:)
        type(Mixture_Short_Potentials_Wrapper), target, intent(in) :: short_potentials
        type(Mixture_Ewald_Wrapper), target, intent(in) :: ewalds
    end subroutine Null_One_Particle_Change_set_candidates

    subroutine Null_One_Particle_Change_set_selector(this)
        class(Null_One_Particle_Change), intent(inout) :: this
    end subroutine Null_One_Particle_Change_set_selector

    subroutine Null_One_Particle_Change_set_observables(this, observables)
        class(Null_One_Particle_Change), intent(inout) :: this
        type(Mixture_Observables_Wrapper), target, intent(in) :: observables
    end subroutine Null_One_Particle_Change_set_observables

    pure integer function Null_One_Particle_Change_get_num_choices(this) result(num_choices)
        class(Null_One_Particle_Change), intent(in) :: this
        num_choices = 0
    end function Null_One_Particle_Change_get_num_choices

    subroutine Null_One_Particle_Change_try(this)
        class(Null_One_Particle_Change), intent(inout) :: this
    end subroutine Null_One_Particle_Change_try

    subroutine Null_One_Particle_Change_define_change(this, new, old, i_actor)
        class(Null_One_Particle_Change), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor
        new%i = 0
        old%i = 0
    end subroutine Null_One_Particle_Change_define_change

    subroutine Null_One_Particle_Change_update_actor(this, i_actor, new, old)
        class(Null_One_Particle_Change), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        new%i = 0
        old%i = 0
    end subroutine Null_One_Particle_Change_update_actor

    subroutine Null_One_Particle_Change_increment_hits(this, i_actor)
        class(Null_One_Particle_Change), intent(inout) :: this
        integer, intent(in) :: i_actor
    end subroutine Null_One_Particle_Change_increment_hits

    subroutine Null_One_Particle_Change_increment_success(this, i_actor)
        class(Null_One_Particle_Change), intent(inout) :: this
        integer, intent(in) :: i_actor
    end subroutine Null_One_Particle_Change_increment_success

!end implementation Null_One_Particle_Change

end module class_one_particle_change
