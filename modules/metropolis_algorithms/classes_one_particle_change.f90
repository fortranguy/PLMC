module classes_one_particle_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_move => visit_move, &
    dipoles_field_visit_rotation => visit_rotation
use classes_tower_sampler, only: Abstract_Tower_Sampler
use module_changes_success, only: Concrete_Changes_Counter
use types_temporary_observables, only: Concrete_Single_Delta_Energies
use types_observables_wrapper, only: Observables_Wrapper
use procedures_metropolis_micro, only: update_energies
use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_One_Particle_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Wrapper), pointer :: dipolar_interactions => null()
        type(Changes_Component_Wrapper), pointer :: change_components(:) => null()
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure :: set_selector => Abstract_set_selector
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure, private :: test_metropolis => Abstract_test_metropolis
        procedure(Abstract_define_change), private, deferred :: define_change
        procedure(Abstract_visit_field), private, deferred :: visit_field
        procedure(Abstract_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_visit_short), private, deferred :: visit_short
        procedure(Abstract_visit_dipolar), private, deferred :: visit_dipolar
        procedure(Abstract_update_actor), private, deferred :: update_actor
        procedure(Abstract_increment_hit), private, nopass, deferred :: increment_hit
        procedure(Abstract_increment_success), private, nopass, deferred :: increment_success
    end type Abstract_One_Particle_Change

    abstract interface

        subroutine Abstract_define_change(this, i_actor, new, old)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(out) :: new, old
        end subroutine Abstract_define_change

        subroutine Abstract_visit_field(this, delta, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            real(DP), intent(out) :: delta
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_field

        subroutine Abstract_visit_walls(this, overlap, delta, i_actor, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, deltas, i_actor, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: deltas(:)
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_short

        subroutine Abstract_visit_dipolar(this, deltas, mixture_delta, i_actor, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            real(DP), intent(out) :: deltas(:)
            real(DP), intent(out) :: mixture_delta
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_actor(this, i_actor, new, old)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Change
            class(Abstract_One_Particle_Change), intent(in) :: this
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_update_actor

        subroutine Abstract_increment_hit(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_increment_hit

        subroutine Abstract_increment_success(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_increment_success

    end interface

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Move
    contains
        procedure, private :: define_change => Move_define_change
        procedure, private :: visit_field => Move_visit_field
        procedure, private :: visit_walls => Move_visit_walls
        procedure, private :: visit_short => Move_visit_short
        procedure, private :: visit_dipolar => Move_visit_dipolar
        procedure, private :: update_actor => Move_update_actor
        procedure, private, nopass :: increment_hit => Move_increment_hit
        procedure, private, nopass :: increment_success => Move_increment_success
    end type Concrete_One_Particle_Move

    type, extends(Abstract_One_Particle_Change), public :: Concrete_One_Particle_Rotation
    contains
        procedure, private :: define_change => Rotation_define_change
        procedure, private :: visit_field => Rotation_visit_field
        procedure, private :: visit_walls => Rotation_visit_walls
        procedure, private :: visit_short => Rotation_visit_short
        procedure, private :: visit_dipolar => Rotation_visit_dipolar
        procedure, private :: update_actor => Rotation_update_actor
        procedure, private, nopass :: increment_hit => Rotation_increment_hit
        procedure, private, nopass :: increment_success => Rotation_increment_success
    end type Concrete_One_Particle_Rotation

    type, extends(Abstract_One_Particle_Change), public :: Null_One_Particle_Change
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set_selector => Null_set_selector
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
        procedure, private :: test_metropolis => Null_test_metropolis
        procedure, private :: define_change => Null_define_change
        procedure, private :: visit_field => Null_visit_field
        procedure, private :: visit_walls => Null_visit_walls
        procedure, private :: visit_short => Null_visit_short
        procedure, private :: visit_dipolar => Null_visit_dipolar
        procedure, private :: update_actor => Null_update_actor
        procedure, private, nopass :: increment_hit => Null_increment_hit
        procedure, private, nopass :: increment_success => Null_increment_success
    end type Null_One_Particle_Change

contains

!implementation Abstract_One_Particle_Change

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions, change_components)
        class(Abstract_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        type(Changes_Component_Wrapper), target, intent(in) :: change_components(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions => dipolar_interactions
        this%change_components => change_components
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_One_Particle_Change), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        this%change_components => null()
        this%dipolar_interactions => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set_selector(this, selector)
        class(Abstract_One_Particle_Change), intent(inout) :: this
        class(Abstract_Tower_Sampler), intent(in) :: selector

        allocate(this%selector, source=selector)
    end subroutine Abstract_set_selector

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_One_Particle_Change), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_One_Particle_Change), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Delta_Energies) :: deltas
        integer :: i_actor

        i_actor = this%selector%get()
        call this%increment_hit(observables%changes_counters(i_actor))
        allocate(deltas%short(size(observables%short_energies)))
        allocate(deltas%dipolar(size(observables%dipolar_energies)))
        call this%test_metropolis(success, deltas, i_actor)
        if (success) then
            observables%field_energies(i_actor) = observables%field_energies(i_actor) + deltas%field
            observables%walls_energies(i_actor) = observables%walls_energies(i_actor) + deltas%walls
            call update_energies(observables%short_energies, deltas%short, i_actor)
            call update_energies(observables%dipolar_energies, deltas%dipolar, i_actor)
            observables%dipolar_mixture_energy = observables%dipolar_mixture_energy + &
                deltas%dipolar_mixture
            call this%increment_success(observables%changes_counters(i_actor))
        end if
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, deltas, i_actor)
        class(Abstract_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Delta_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_actor

        type(Concrete_Temporary_Particle) :: new, old
        real(DP) :: energy_delta
        logical :: overlap
        real(DP) :: rand

        call this%define_change(i_actor, new, old)

        success = .false.
        call this%visit_field(deltas%field, new, old)
        call this%visit_walls(overlap, deltas%walls, i_actor, new, old)
        if (overlap) return
        call this%visit_short(overlap, deltas%short, i_actor, new, old)
        if (overlap) return
        call this%visit_dipolar(deltas%dipolar, deltas%dipolar_mixture, i_actor, new, old)

        energy_delta = deltas%field + deltas%walls + sum(deltas%short + deltas%dipolar) + &
            deltas%dipolar_mixture
        call random_number(rand)
        if (rand < exp(-energy_delta/this%environment%temperature%get())) then
            call this%update_actor(i_actor, new, old)
            success = .true.
        end if
    end subroutine Abstract_test_metropolis

!end implementation Abstract_One_Particle_Change

!implementation Concrete_One_Particle_Move

    subroutine Move_define_change(this, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(out) :: new, old

        old%i = random_integer(this%mixture%components(i_actor)%positions%get_num())
        old%position = this%mixture%components(i_actor)%positions%get(old%i)
        old%orientation = this%mixture%components(i_actor)%orientations%get(old%i)
        old%dipolar_moment = this%mixture%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = this%change_components(i_actor)%moved_positions%get(new%i)
        new%orientation = old%orientation
        new%dipolar_moment = old%dipolar_moment
    end subroutine Move_define_change

    subroutine Move_visit_field(this, delta, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        real(DP), intent(out) :: delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta = dipoles_field_visit_move(this%environment%external_field, new%position, old)
    end subroutine Move_visit_field

    subroutine Move_visit_walls(this, overlap, delta, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: energy_new, energy_old

        call this%environment%walls%visit(overlap, energy_new, new%position, this%&
            short_interactions%wall_pairs(i_actor)%potential)
        if (overlap) return
        call this%environment%walls%visit(overlap, energy_old, old%position, this%&
            short_interactions%wall_pairs(i_actor)%potential)
        delta = energy_new - energy_old
    end subroutine Move_visit_walls

    subroutine Move_visit_short(this, overlap, deltas, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(deltas)) :: energies_new, energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%short_interactions%visitable_cells, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%short_interactions%visitable_cells(i_component, i_actor)%visit(overlap, &
                energies_new(i_component), new, i_exclude)
            if (overlap) return
        end do
        do i_component = 1, size(this%short_interactions%visitable_cells, 1)
            i_exclude = merge(old%i, 0, i_component == i_actor)
            call this%short_interactions%visitable_cells(i_component, i_actor)%visit(overlap, &
                energies_old(i_component), old, i_exclude)
        end do
        deltas = energies_new - energies_old
    end subroutine Move_visit_short

    subroutine Move_visit_dipolar(this, deltas, mixture_delta, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        real(DP), intent(out) :: mixture_delta
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(deltas)) :: real_energies_new, real_energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%dipolar_interactions%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%dipolar_interactions%real_components(i_component, i_actor)%component%&
                visit(real_energies_new(i_component), new, i_exclude)
            call this%dipolar_interactions%real_components(i_component, i_actor)%component%&
                visit(real_energies_old(i_component), old, i_exclude)
        end do
        mixture_delta = this%dipolar_interactions%reci_visitor%visit_move(i_actor, new%position, &
            old) - this%dipolar_interactions%dlc_visitor%visit_move(i_actor, new%position, old)
        deltas = real_energies_new - real_energies_old
    end subroutine Move_visit_dipolar

    subroutine Move_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Move), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        integer :: i_component

        call this%mixture%components(i_actor)%positions%set(new%i, new%position)
        do i_component = 1, size(this%short_interactions%visitable_cells, 1)
            call this%short_interactions%visitable_cells(i_actor, i_component)%move(new%position, &
                old)
        end do
        call this%dipolar_interactions%reci_structure%update_move(i_actor, new%position, old)
        call this%dipolar_interactions%dlc_structures%update_move(i_actor, new%position, old)
    end subroutine Move_update_actor

    subroutine Move_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%move%num_hits = changes_counters%move%num_hits + 1
    end subroutine Move_increment_hit

    subroutine Move_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%move%num_success = changes_counters%move%num_success + 1
    end subroutine Move_increment_success

!end implementation Concrete_One_Particle_Move

!implementation Concrete_One_Particle_Rotation

    subroutine Rotation_define_change(this, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(out) :: new, old

        old%i = random_integer(this%mixture%components(i_actor)%orientations%get_num())
        old%position = this%mixture%components(i_actor)%positions%get(old%i)
        old%orientation = this%mixture%components(i_actor)%orientations%get(old%i)
        old%dipolar_moment = this%mixture%components(i_actor)%dipolar_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%change_components(i_actor)%rotated_orientations%get(new%i)
        new%dipolar_moment = this%mixture%components(i_actor)%dipolar_moments%get_norm() * new%&
            orientation
    end subroutine Rotation_define_change

    subroutine Rotation_visit_field(this, delta, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta = dipoles_field_visit_rotation(this%environment%external_field, new%dipolar_moment, &
            old)
    end subroutine Rotation_visit_field

    subroutine Rotation_visit_walls(this, overlap, delta, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        delta = 0._DP
    end subroutine Rotation_visit_walls

    subroutine Rotation_visit_short(this, overlap, deltas, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        deltas = 0._DP
    end subroutine Rotation_visit_short

    subroutine Rotation_visit_dipolar(this, deltas, mixture_delta, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        real(DP), intent(out) :: mixture_delta
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(deltas)) :: real_energies_new, real_energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%dipolar_interactions%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%dipolar_interactions%real_components(i_component, i_actor)%component%&
                visit(real_energies_new(i_component), new, i_exclude)
            call this%dipolar_interactions%real_components(i_component, i_actor)%component%&
                visit(real_energies_old(i_component), old, i_exclude)
        end do
        deltas = real_energies_new - real_energies_old
        mixture_delta = this%dipolar_interactions%reci_visitor%visit_rotation(i_actor, new%&
            dipolar_moment, old) + this%dipolar_interactions%surf_mixture%visit_rotation(i_actor, &
            new%dipolar_moment, old%dipolar_moment) - this%dipolar_interactions%dlc_visitor%&
            visit_rotation(i_actor, new%dipolar_moment, old)
    end subroutine Rotation_visit_dipolar

    subroutine Rotation_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        call this%mixture%components(i_actor)%orientations%set(new%i, new%orientation)
        call this%mixture%total_moment%remove(i_actor, old%dipolar_moment)
        call this%mixture%total_moment%add(i_actor, new%dipolar_moment)
        call this%dipolar_interactions%reci_structure%update_rotation(i_actor, new%dipolar_moment, &
            old)
        call this%dipolar_interactions%dlc_structures%update_rotation(i_actor, new%dipolar_moment, &
            old)
    end subroutine Rotation_update_actor

    subroutine Rotation_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_hits = changes_counters%rotation%num_hits + 1
    end subroutine Rotation_increment_hit

    subroutine Rotation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_success = changes_counters%rotation%num_success + 1
    end subroutine Rotation_increment_success

!end implementation Concrete_One_Particle_Rotation

!implementation Null_One_Particle_Change

    subroutine Null_construct(this, environment, mixture, short_interactions, dipolar_interactions,&
        change_components)
        class(Null_One_Particle_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        type(Changes_Component_Wrapper), target, intent(in) :: change_components(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_One_Particle_Change), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_selector(this, selector)
        class(Null_One_Particle_Change), intent(inout) :: this
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_set_selector

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_One_Particle_Change), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_One_Particle_Change), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

    subroutine Null_test_metropolis(this, success, deltas, i_actor)
        class(Null_One_Particle_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Delta_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_actor
        success = .false.
        deltas%field = 0._DP; deltas%walls = 0._DP; deltas%short = 0._DP; deltas%dipolar = 0._DP
        deltas%dipolar_mixture = 0._DP
    end subroutine Null_test_metropolis

    subroutine Null_define_change(this, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(out) :: new, old
         new%i = 0; new%position = 0._DP; new%orientation = 0._DP; new%dipolar_moment = 0._DP
         old%i = 0; old%position = 0._DP; old%orientation = 0._DP; old%dipolar_moment = 0._DP
     end subroutine Null_define_change

     subroutine Null_visit_field(this, delta, new, old)
        class(Null_One_Particle_Change), intent(in) :: this
        real(DP), intent(out) :: delta
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        delta = 0._DP
    end subroutine Null_visit_field

     subroutine Null_visit_walls(this, overlap, delta, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: delta
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         overlap = .false.
         delta = 0._DP
     end subroutine Null_visit_walls

     subroutine Null_visit_short(this, overlap, deltas, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: deltas(:)
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         overlap = .false.
         deltas = 0._DP
     end subroutine Null_visit_short

     subroutine Null_visit_dipolar(this, deltas, mixture_delta, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         real(DP), intent(out) :: deltas(:)
         real(DP), intent(out) :: mixture_delta
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         integer, intent(in) :: i_actor
         deltas = 0._DP
         mixture_delta = 0._DP
     end subroutine Null_visit_dipolar

     subroutine Null_update_actor(this, i_actor, new, old)
         class(Null_One_Particle_Change), intent(in) :: this
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
     end subroutine Null_update_actor

    subroutine Null_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_increment_hit

    subroutine Null_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_increment_success

!end implementation Null_One_Particle_Change

end module classes_one_particle_change
