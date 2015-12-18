module class_one_particle_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler
use module_changes_success, only: Concrete_Changes_Counter
use types_observables_wrapper, only: Observables_Wrapper
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_One_Particle_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Changes_Component_Wrapper), pointer :: changes(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Long_Interactions_Wrapper), pointer :: long_interactions => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_One_Particle_Change_construct
        procedure :: destroy => Abstract_One_Particle_Change_destroy
        procedure :: try => Abstract_One_Particle_Change_try
        procedure :: set_candidates => Abstract_One_Particle_Change_set_candidates
        procedure :: get_num_choices => Abstract_One_Particle_Change_get_num_choices
        procedure, private :: test_metropolis => Abstract_One_Particle_Change_test_metropolis
        procedure(Abstract_One_Particle_Change_define_change), private, deferred :: define_change
        procedure(Abstract_One_Particle_Change_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_One_Particle_Change_visit_short), private, deferred :: visit_short
        procedure(Abstract_One_Particle_Change_visit_long), private, deferred :: visit_long
        procedure(Abstract_One_Particle_Change_update_actor), private, deferred :: update_actor
        procedure(Abstract_One_Particle_Change_increment_hits), private, nopass, deferred :: &
            increment_hits
        procedure(Abstract_One_Particle_Change_increment_success), private, nopass, deferred :: &
            increment_success
    end type Abstract_One_Particle_Change

    abstract interface

        subroutine Abstract_One_Particle_Change_define_change(this, new, old, i_actor)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            type(Concrete_Temporary_Particle), intent(out) :: new, old
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_define_change

        subroutine Abstract_One_Particle_Change_visit_walls(this, overlap, walls_delta, new, old, &
            i_actor)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: walls_delta
            type(Concrete_Temporary_Particle), intent(in) :: new, old
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_visit_walls

        subroutine Abstract_One_Particle_Change_visit_short(this, overlap, short_deltas, new, old, &
            i_actor)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            logical, intent(out) :: overlap
            type(Concrete_Temporary_Particle), intent(in) :: new, old
            real(DP), intent(out) :: short_deltas(:)
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_visit_short

        subroutine Abstract_One_Particle_Change_visit_long(this, long_deltas, reci_delta, new, &
            old, i_actor)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            real(DP), intent(out) :: long_deltas(:)
            real(DP), intent(out) :: reci_delta
            type(Concrete_Temporary_Particle), intent(in) :: new, old
            integer, intent(in) :: i_actor
        end subroutine Abstract_One_Particle_Change_visit_long

        subroutine Abstract_One_Particle_Change_update_actor(this, i_actor, new, old)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_One_Particle_Change_update_actor

        subroutine Abstract_One_Particle_Change_increment_hits(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_One_Particle_Change_increment_hits

        subroutine Abstract_One_Particle_Change_increment_success(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_One_Particle_Change_increment_success

    end interface

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Move
    contains
        procedure, private :: define_change => Concrete_One_Particle_Move_define_change
        procedure, private :: visit_walls => Concrete_One_Particle_Move_visit_walls
        procedure, private :: visit_short => Concrete_One_Particle_Move_visit_short
        procedure, private :: visit_long => Concrete_One_Particle_Move_visit_long
        procedure, private :: update_actor => Concrete_One_Particle_Move_update_actor
        procedure, private, nopass :: increment_hits => Concrete_One_Particle_Move_increment_hits
        procedure, private, nopass :: increment_success => &
            Concrete_One_Particle_Move_increment_success
    end type Concrete_One_Particle_Move

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Rotation
    contains
        procedure, private :: define_change => Concrete_One_Particle_Rotation_define_change
        procedure, private :: visit_walls => Concrete_One_Particle_Rotation_visit_walls
        procedure, private :: visit_short => Concrete_One_Particle_Rotation_visit_short
        procedure, private :: visit_long => Concrete_One_Particle_Rotation_visit_long
        procedure, private :: update_actor => Concrete_One_Particle_Rotation_update_actor
        procedure, private, nopass :: increment_hits => &
            Concrete_One_Particle_Rotation_increment_hits
        procedure, private, nopass :: increment_success => &
            Concrete_One_Particle_Rotation_increment_success
    end type Concrete_One_Particle_Rotation

    type, extends(Abstract_One_Particle_Change), public :: Null_One_Particle_Change
    contains
        procedure :: construct => Null_One_Particle_Change_construct
        procedure :: destroy => Null_One_Particle_Change_destroy
        procedure :: set_candidates => Null_One_Particle_Change_set_candidates
        procedure :: get_num_choices => Null_One_Particle_Change_get_num_choices
        procedure :: try => Null_One_Particle_Change_try
        procedure, private :: test_metropolis => Null_One_Particle_Change_test_metropolis
        procedure, private :: define_change => Null_One_Particle_Change_define_change
        procedure, private :: visit_walls => Null_One_Particle_Change_visit_walls
        procedure, private :: visit_short => Null_One_Particle_Change_visit_short
        procedure, private :: visit_long => Null_One_Particle_Change_visit_long
        procedure, private :: update_actor => Null_One_Particle_Change_update_actor
        procedure, private, nopass :: increment_hits => Null_One_Particle_Change_increment_hits
        procedure, private, nopass :: increment_success => &
            Null_One_Particle_Change_increment_success
    end type Null_One_Particle_Change

contains

!implementation Abstract_One_Particle_Change

    subroutine Abstract_One_Particle_Change_construct(this, environment, changes, selector)
        class(Abstract_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Changes_Component_Wrapper), target, intent(in) :: changes(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%environment => environment
        this%changes => changes
        allocate(this%selector, source=selector)
    end subroutine Abstract_One_Particle_Change_construct

    subroutine Abstract_One_Particle_Change_destroy(this)
        class(Abstract_One_Particle_Change), intent(inout) :: this

        this%long_interactions => null()
        this%short_interactions => null()
        this%components => null()
        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        this%changes => null()
        this%environment => null()
    end subroutine Abstract_One_Particle_Change_destroy

    subroutine Abstract_One_Particle_Change_set_candidates(this, components, short_interactions, &
        long_interactions)
        class(Abstract_One_Particle_Change), intent(inout) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), target, intent(in) :: long_interactions

        this%components => components
        this%short_interactions => short_interactions
        this%long_interactions => long_interactions
    end subroutine Abstract_One_Particle_Change_set_candidates

    pure integer function Abstract_One_Particle_Change_get_num_choices(this) result(num_choices)
        class(Abstract_One_Particle_Change), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_One_Particle_Change_get_num_choices

    subroutine Abstract_One_Particle_Change_try(this, observables)
        class(Abstract_One_Particle_Change), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        integer :: i_actor, i_component
        logical :: success
        real(DP) :: walls_delta, field_energy
        real(DP) :: short_deltas(size(observables%short_energies)), long_deltas(size(observables%&
            long_energies_wo_reci))
        real(DP) :: reci_delta
        integer :: j_observable, i_observable

        i_actor = this%selector%get()
        call this%increment_hits(observables%changes_counters(i_actor))
        call this%test_metropolis(success, walls_delta, short_deltas, long_deltas, reci_delta, &
            i_actor)
        if (success) then
            do i_component = 1, size(observables%short_energies)
                j_observable = maxval([i_actor, i_component])
                i_observable = minval([i_actor, i_component])
                observables%short_energies(j_observable)%with_components(i_observable) = &
                    observables%short_energies(j_observable)%with_components(i_observable) + &
                    short_deltas(i_component)
            end do
            do i_component = 1, size(observables%long_energies_wo_reci)
                j_observable = maxval([i_actor, i_component])
                i_observable = minval([i_actor, i_component])
                observables%long_energies_wo_reci(j_observable)%with_components(i_observable) = &
                    observables%long_energies_wo_reci(j_observable)%with_components(i_observable) +&
                    long_deltas(i_component)
            end do
            observables%reci_energy = observables%reci_energy + reci_delta
            call this%increment_success(observables%changes_counters(i_actor))
        end if
    end subroutine Abstract_One_Particle_Change_try

    subroutine Abstract_One_Particle_Change_test_metropolis(this, success, walls_delta, &
        short_deltas, long_deltas, reci_delta, i_actor)
        class(Abstract_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: success
        real(DP), intent(out) :: walls_delta
        real(DP), intent(out) :: short_deltas(:), long_deltas(:)
        real(DP), intent(out) :: reci_delta
        integer, intent(in) :: i_actor

        type(Concrete_Temporary_Particle) :: new, old
        real(DP) :: energy_delta
        logical :: overlap
        real(DP) :: rand

        call this%define_change(new, old, i_actor)

        success = .false.
        call this%visit_walls(overlap, walls_delta, new, old, i_actor)
        if (overlap) return
        call this%visit_short(overlap, short_deltas, new, old, i_actor)
        if (overlap) return
        call this%visit_long(long_deltas, reci_delta, new, old, i_actor)

        energy_delta = walls_delta + sum(short_deltas + long_deltas) + reci_delta
        call random_number(rand)
        if (rand < exp(-energy_delta/this%environment%temperature%get())) then
            call this%update_actor(i_actor, new, old)
            success = .true.
        end if
    end subroutine Abstract_One_Particle_Change_test_metropolis

!end implementation Abstract_One_Particle_Change

!implementation Concrete_One_Particle_Move

    subroutine Concrete_One_Particle_Move_define_change(this, new, old, i_actor)
        class(Concrete_One_Particle_Move), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        old%i = random_integer(this%components(i_actor)%positions%get_num())
        old%position = this%components(i_actor)%positions%get(old%i)
        old%orientation = this%components(i_actor)%orientations%get(old%i)
        old%dipolar_moment = this%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = this%changes(i_actor)%moved_positions%get(new%i)
        new%orientation = old%orientation
        new%dipolar_moment = old%dipolar_moment
    end subroutine Concrete_One_Particle_Move_define_change

    subroutine Concrete_One_Particle_Move_visit_walls(this, overlap, walls_delta, new, old, i_actor)
        class(Concrete_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: walls_delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor

        real(DP) :: walls_new, walls_old

        call this%environment%walls_potential%visit(overlap, walls_new, new%position, &
            this%short_interactions%wall_pairs(i_actor)%pair_potential)
        if (overlap) return
        call this%environment%walls_potential%visit(overlap, walls_old, old%position, &
            this%short_interactions%wall_pairs(i_actor)%pair_potential)
        walls_delta = walls_new - walls_old
    end subroutine Concrete_One_Particle_Move_visit_walls

    subroutine Concrete_One_Particle_Move_visit_short(this, overlap, short_deltas, new, old, &
        i_actor)
        class(Concrete_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: short_deltas(:)
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor

        real(DP), dimension(size(short_deltas)) :: short_new, short_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%short_interactions%components_cells, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%short_interactions%components_cells(i_component, i_actor)%visit(overlap, &
                short_new(i_component), new, i_exclude)
            if (overlap) return
        end do
        do i_component = 1, size(this%short_interactions%components_cells, 1)
            i_exclude = merge(old%i, 0, i_component == i_actor)
            call this%short_interactions%components_cells(i_component, i_actor)%visit(overlap, &
                short_old(i_component), old, i_exclude)
        end do
        short_deltas = short_new - short_old
    end subroutine Concrete_One_Particle_Move_visit_short

    subroutine Concrete_One_Particle_Move_visit_long(this, long_deltas, reci_delta, new, old, &
        i_actor)
        class(Concrete_One_Particle_Move), intent(in) :: this
        real(DP), intent(out) :: long_deltas(:)
        real(DP), intent(out) :: reci_delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor

        real(DP), dimension(size(long_deltas)) :: long_new_real, long_old_real
        integer :: i_component, i_exclude

        do i_component = 1, size(this%long_interactions%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%long_interactions%real_components(i_component, i_actor)%real_component%&
                visit(long_new_real(i_component), new, i_exclude)
            call this%long_interactions%real_components(i_component, i_actor)%real_component%&
                visit(long_old_real(i_component), old, i_exclude)
        end do
        reci_delta = this%long_interactions%reci_visitor%visit_move(i_actor, new, old)
        long_deltas = long_new_real-long_old_real
    end subroutine Concrete_One_Particle_Move_visit_long

    subroutine Concrete_One_Particle_Move_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        integer :: i_component

        call this%components(i_actor)%positions%set(new%i, new%position)
        do i_component = 1, size(this%short_interactions%components_cells, 1)
            call this%short_interactions%components_cells(i_actor, i_component)%move(old, new)
        end do
        call this%long_interactions%reci_structure%update_move(i_actor, new, old)
    end subroutine Concrete_One_Particle_Move_update_actor

    subroutine Concrete_One_Particle_Move_increment_hits(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%move%num_hits = changes_counters%move%num_hits + 1
    end subroutine Concrete_One_Particle_Move_increment_hits

    subroutine Concrete_One_Particle_Move_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%move%num_success = changes_counters%move%num_success + 1
    end subroutine Concrete_One_Particle_Move_increment_success

!end implementation Concrete_One_Particle_Move

!implementation Concrete_One_Particle_Rotation

    subroutine Concrete_One_Particle_Rotation_define_change(this, new, old, i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        old%i = random_integer(this%components(i_actor)%orientations%get_num())
        old%position = this%components(i_actor)%positions%get(old%i)
        old%orientation = this%components(i_actor)%orientations%get(old%i)
        old%dipolar_moment = this%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%changes(i_actor)%rotated_orientations%get(new%i)
        new%dipolar_moment = this%components(i_actor)%dipolar_moments%get_norm() * new%orientation
    end subroutine Concrete_One_Particle_Rotation_define_change

    subroutine Concrete_One_Particle_Rotation_visit_walls(this, overlap, walls_delta, new, &
        old, i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: walls_delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor
        overlap = .false.
        walls_delta = 0._DP
    end subroutine Concrete_One_Particle_Rotation_visit_walls

    subroutine Concrete_One_Particle_Rotation_visit_short(this, overlap, short_deltas, new, &
        old, i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: short_deltas(:)
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor
        overlap = .false.
        short_deltas = 0._DP
    end subroutine Concrete_One_Particle_Rotation_visit_short

    subroutine Concrete_One_Particle_Rotation_visit_long(this, long_deltas, reci_delta, new, old, &
        i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: long_deltas(:)
        real(DP), intent(out) :: reci_delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        integer, intent(in) :: i_actor

        real(DP), dimension(size(long_deltas)) :: long_new_real, long_old_real
        real(DP) :: long_delta_self
        integer :: i_component, i_exclude

        do i_component = 1, size(this%long_interactions%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%long_interactions%real_components(i_component, i_actor)%real_component%&
                visit(long_new_real(i_component), new, i_exclude)
            call this%long_interactions%real_components(i_component, i_actor)%real_component%&
                visit(long_old_real(i_component), old, i_exclude)
        end do
        reci_delta = this%long_interactions%reci_visitor%visit_rotation(i_actor, new, old)
        long_delta_self = this%long_interactions%self_components(i_actor)%self%&
            meet(new%dipolar_moment) - this%long_interactions%self_components(i_actor)%self%&
            meet(old%dipolar_moment)
        long_deltas = long_new_real-long_old_real - long_delta_self
    end subroutine Concrete_One_Particle_Rotation_visit_long

    subroutine Concrete_One_Particle_Rotation_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        call this%components(i_actor)%orientations%set(new%i, new%orientation)
        call this%long_interactions%reci_structure%update_rotation(i_actor, new, old)
    end subroutine Concrete_One_Particle_Rotation_update_actor

    subroutine Concrete_One_Particle_Rotation_increment_hits(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_hits = changes_counters%rotation%num_hits + 1
    end subroutine Concrete_One_Particle_Rotation_increment_hits

    subroutine Concrete_One_Particle_Rotation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_success = changes_counters%rotation%num_success + 1
    end subroutine Concrete_One_Particle_Rotation_increment_success

!end implementation Concrete_One_Particle_Rotation

!implementation Null_One_Particle_Change

    subroutine Null_One_Particle_Change_construct(this, environment, changes, selector)
        class(Null_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Changes_Component_Wrapper), target, intent(in) :: changes(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_One_Particle_Change_construct

    subroutine Null_One_Particle_Change_destroy(this)
        class(Null_One_Particle_Change), intent(inout) :: this
    end subroutine Null_One_Particle_Change_destroy

    subroutine Null_One_Particle_Change_set_candidates(this, components, short_interactions, &
        long_interactions)
        class(Null_One_Particle_Change), intent(inout) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), target, intent(in) :: long_interactions
    end subroutine Null_One_Particle_Change_set_candidates

    pure integer function Null_One_Particle_Change_get_num_choices(this) result(num_choices)
        class(Null_One_Particle_Change), intent(in) :: this
        num_choices = 0
    end function Null_One_Particle_Change_get_num_choices

    subroutine Null_One_Particle_Change_try(this, observables)
        class(Null_One_Particle_Change), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_One_Particle_Change_try

    subroutine Null_One_Particle_Change_test_metropolis(this, success, walls_delta, &
        short_deltas, long_deltas, reci_delta, i_actor)
        class(Null_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: success
        real(DP), intent(out) :: walls_delta
        real(DP), intent(out) :: short_deltas(:), long_deltas(:)
        real(DP), intent(out) :: reci_delta
        integer, intent(in) :: i_actor
        success = .false.
        walls_delta = 0._DP; short_deltas = 0._DP; long_deltas = 0._DP; reci_delta = 0._DP
    end subroutine Null_One_Particle_Change_test_metropolis

    subroutine Null_One_Particle_Change_define_change(this, new, old, i_actor)
         class(Null_One_Particle_Change), intent(in) :: this
         type(Concrete_Temporary_Particle), intent(out) :: new, old
         integer, intent(in) :: i_actor
         new%i = 0; new%position = 0._DP; new%orientation = 0._DP; new%dipolar_moment = 0._DP
         old%i = 0; old%position = 0._DP; old%orientation = 0._DP; old%dipolar_moment = 0._DP
     end subroutine Null_One_Particle_Change_define_change

     subroutine Null_One_Particle_Change_visit_walls(this, overlap, walls_delta, new, old, i_actor)
         class(Null_One_Particle_Change), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: walls_delta
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         integer, intent(in) :: i_actor
         overlap = .false.
         walls_delta = 0._DP
     end subroutine Null_One_Particle_Change_visit_walls

     subroutine Null_One_Particle_Change_visit_short(this, overlap, short_deltas, new, old, i_actor)
         class(Null_One_Particle_Change), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: short_deltas(:)
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         integer, intent(in) :: i_actor
         overlap = .false.
         short_deltas = 0._DP
     end subroutine Null_One_Particle_Change_visit_short

     subroutine Null_One_Particle_Change_visit_long(this, long_deltas, reci_delta, new, old, &
        i_actor)
         class(Null_One_Particle_Change), intent(in) :: this
         real(DP), intent(out) :: long_deltas(:)
         real(DP), intent(out) :: reci_delta
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         integer, intent(in) :: i_actor
         long_deltas = 0._DP
         reci_delta = 0._DP
     end subroutine Null_One_Particle_Change_visit_long

     subroutine Null_One_Particle_Change_update_actor(this, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
     end subroutine Null_One_Particle_Change_update_actor

    subroutine Null_One_Particle_Change_increment_hits(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_One_Particle_Change_increment_hits

    subroutine Null_One_Particle_Change_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_One_Particle_Change_increment_success

!end implementation Null_One_Particle_Change

end module class_one_particle_change
