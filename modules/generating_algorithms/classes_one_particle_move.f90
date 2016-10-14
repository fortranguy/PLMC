module classes_one_particle_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_different
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipoles_field_interaction, only: &
    dipoles_field_visit_translation => visit_translation, &
    dipoles_field_visit_rotation => visit_rotation
use types_changes_component_wrapper, only: Changes_Component_Wrapper
use module_changes_success, only: Concrete_Changes_Counter
use types_observables_energies, only: Concrete_Single_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Abstract_One_Particle_Move
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static => null()
        type(Changes_Component_Wrapper), pointer :: changes_components(:) => null()
        logical, allocatable :: can_move(:)
        class(Abstract_Tower_Sampler), allocatable :: selector
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure :: reset_selector => Abstract_reset_selector
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure(Abstract_define_change), private, deferred :: define_change
        procedure(Abstract_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_visit_short), private, deferred :: visit_short
        procedure(Abstract_visit_field), private, deferred :: visit_field
        procedure(Abstract_visit_dipolar), private, deferred :: visit_dipolar
        procedure(Abstract_update_actor), private, deferred :: update_actor
        procedure(Abstract_increment_hit), private, nopass, deferred :: increment_hit
        procedure(Abstract_increment_success), private, nopass, deferred :: increment_success
    end type Abstract_One_Particle_Move

    abstract interface

        subroutine Abstract_define_change(this, abort, new, old, i_actor)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
            logical, intent(out) :: abort
            type(Concrete_Temporary_Particle), intent(out) :: new, old
            integer, intent(in) :: i_actor
        end subroutine Abstract_define_change

        subroutine Abstract_visit_walls(this, overlap, delta_energy, i_actor, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, delta_energies, i_actor, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_short

        subroutine Abstract_visit_field(this, delta_energy, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
            real(DP), intent(out) :: delta_energy
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_field

        subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, i_actor, new, &
            old)
        import :: DP, Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:)
            real(DP), intent(out) :: delta_shared_energy
            integer, intent(in) :: i_actor
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_actor(this, i_actor, new, old)
        import :: Concrete_Temporary_Particle, Abstract_One_Particle_Move
            class(Abstract_One_Particle_Move), intent(in) :: this
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

    type, extends(Abstract_One_Particle_Move), public :: Concrete_One_Particle_Translation
    contains
        procedure, private :: define_change => Translation_define_change
        procedure, private :: visit_walls => Translation_visit_walls
        procedure, private :: visit_short => Translation_visit_short
        procedure, private :: visit_field => Translation_visit_field
        procedure, private :: visit_dipolar => Translation_visit_dipolar
        procedure, private :: update_actor => Translation_update_actor
        procedure, private, nopass :: increment_hit => Translation_increment_hit
        procedure, private, nopass :: increment_success => Translation_increment_success
    end type Concrete_One_Particle_Translation

    type, extends(Abstract_One_Particle_Move), public :: Concrete_One_Particle_Rotation
    contains
        procedure, private :: define_change => Rotation_define_change
        procedure, private :: visit_walls => Rotation_visit_walls
        procedure, private :: visit_short => Rotation_visit_short
        procedure, private :: visit_field => Rotation_visit_field
        procedure, private :: visit_dipolar => Rotation_visit_dipolar
        procedure, private :: update_actor => Rotation_update_actor
        procedure, private, nopass :: increment_hit => Rotation_increment_hit
        procedure, private, nopass :: increment_success => Rotation_increment_success
    end type Concrete_One_Particle_Rotation

    type, extends(Abstract_One_Particle_Move), public :: Null_One_Particle_Move
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset_selector => Null_reset_selector
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
        procedure, private :: metropolis_algorithm => Null_metropolis_algorithm
        procedure, private :: define_change => Null_define_change
        procedure, private :: visit_walls => Null_visit_walls
        procedure, private :: visit_short => Null_visit_short
        procedure, private :: visit_field => Null_visit_field
        procedure, private :: visit_dipolar => Null_visit_dipolar
        procedure, private :: update_actor => Null_update_actor
        procedure, private, nopass :: increment_hit => Null_increment_hit
        procedure, private, nopass :: increment_success => Null_increment_success
    end type Null_One_Particle_Move

contains

!implementation Abstract_One_Particle_Move

    !> @note this%selector construction is delayed in
    !> [[classes_one_particle_move:Abstract_reset_selector]]
    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes_components, can_move, &
        selector)
        class(Abstract_One_Particle_Move), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static
        type(Changes_Component_Wrapper), target, intent(in) :: changes_components(:)
        logical, intent(in) :: can_move(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%changes_components => changes_components
        allocate(this%can_move, source=can_move)
        allocate(this%selector, source=selector)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_One_Particle_Move), intent(inout) :: this

        call tower_sampler_destroy(this%selector)
        if (allocated(this%can_move)) deallocate(this%can_move)
        this%changes_components => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selector(this)
        class(Abstract_One_Particle_Move), intent(inout) :: this

        integer :: nums_candidates(size(this%can_move)), i_component

        do i_component = 1, size(nums_candidates)
            nums_candidates(i_component) = merge(this%mixture%components(i_component)%&
                average_num_particles%get(), 0, this%can_move(i_component))
        end do
        call this%selector%reset(nums_candidates)
    end subroutine Abstract_reset_selector

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_One_Particle_Move), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_One_Particle_Move), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Energies) :: deltas
        integer :: i_actor

        i_actor = this%selector%get()
        call this%increment_hit(observables%changes_counters(i_actor))
        allocate(deltas%short_energies(size(observables%energies%short_energies)))
        allocate(deltas%dipolar_energies(size(observables%energies%dipolar_energies)))
        call this%metropolis_algorithm(success, deltas, i_actor)
        if (success) then
            call observables_energies_set(observables%energies, deltas, i_actor)
            call this%increment_success(observables%changes_counters(i_actor))
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, i_actor)
        class(Abstract_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_actor

        real(DP) :: delta_energy
        type(Concrete_Temporary_Particle) :: new, old
        logical :: abort, overlap

        success = .false.
        call this%define_change(abort, new, old, i_actor)
        if (abort) return

        call this%visit_walls(overlap, deltas%walls_energy, i_actor, new, old)
        if (overlap) return
        call this%visit_short(overlap, deltas%short_energies, i_actor, new, old)
        if (overlap) return
        call this%visit_field(deltas%field_energy, new, old)
        call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, i_actor, &
            new, old)

        delta_energy = deltas%walls_energy + deltas%field_energy + &
            sum(deltas%short_energies + deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(min(1._DP, &
            exp(-delta_energy/this%environment%temperature%get())))
        if (success) call this%update_actor(i_actor, new, old)
    end subroutine Abstract_metropolis_algorithm

!end implementation Abstract_One_Particle_Move

!implementation Concrete_One_Particle_Translation

    subroutine Translation_define_change(this, abort, new, old, i_actor)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        if (this%mixture%components(i_actor)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        old%i = random_integer(this%mixture%components(i_actor)%num_particles%get())
        old%position = this%mixture%components(i_actor)%positions%get(old%i)
        old%orientation = this%mixture%components(i_actor)%orientations%get(old%i)
        old%dipole_moment = this%mixture%components(i_actor)%dipole_moments%get(old%i)
        new%i = old%i
        new%position = this%changes_components(i_actor)%translated_positions%get(new%i)
        new%orientation = old%orientation
        new%dipole_moment = old%dipole_moment
    end subroutine Translation_define_change

    subroutine Translation_visit_walls(this, overlap, delta_energy, i_actor, new, old)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: energy_new, energy_old

        call this%environment%visitable_walls%visit(overlap, energy_new, new%position, this%&
            short_interactions%wall_pairs(i_actor)%potential)
        if (overlap) return
        call this%environment%visitable_walls%visit(overlap, energy_old, old%position, this%&
            short_interactions%wall_pairs(i_actor)%potential)
        delta_energy = energy_new - energy_old
    end subroutine Translation_visit_walls

    subroutine Translation_visit_short(this, overlap, delta_energies, i_actor, new, old)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: energies_new, energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%short_interactions%visitable_cells, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%short_interactions%visitable_cells(i_component, i_actor)%&
                visit_energy(overlap, energies_new(i_component), new, visit_different, i_exclude)
            if (overlap) return
        end do
        do i_component = 1, size(this%short_interactions%visitable_cells, 1)
            i_exclude = merge(old%i, 0, i_component == i_actor)
            call this%short_interactions%visitable_cells(i_component, i_actor)%&
                visit_energy(overlap, energies_old(i_component), old, visit_different, i_exclude)
        end do
        delta_energies = energies_new - energies_old
    end subroutine Translation_visit_short

    subroutine Translation_visit_field(this, delta_energy, new, old)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta_energy = dipoles_field_visit_translation(this%environment%external_field, new%&
            position, old)
    end subroutine Translation_visit_field

    subroutine Translation_visit_dipolar(this, delta_energies, delta_shared_energy, i_actor, new, &
        old)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: real_energies_new, real_energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%dipolar_interactions_dynamic%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%dipolar_interactions_dynamic%real_components(i_component, i_actor)%component%&
                visit(real_energies_new(i_component), new, visit_different, i_exclude)
            call this%dipolar_interactions_dynamic%real_components(i_component, i_actor)%component%&
                visit(real_energies_old(i_component), old, visit_different, i_exclude)
        end do
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitor%&
                visit_translation(i_actor, new%position, old) - &
            this%dipolar_interactions_dynamic%dlc_visitor%&
                visit_translation(i_actor, new%position, old)
        delta_energies = real_energies_new - real_energies_old
    end subroutine Translation_visit_dipolar

    subroutine Translation_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Translation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        integer :: i_component

        call this%mixture%components(i_actor)%positions%set(new%i, new%position)
        do i_component = 1, size(this%short_interactions%visitable_cells, 2)
            call this%short_interactions%visitable_cells(i_actor, i_component)%translate(new%&
                position, old)
        end do
        call this%dipolar_interactions_static%reci_structure%&
            update_translation(i_actor, new%position, old)
        call this%dipolar_interactions_static%dlc_structures%&
            update_translation(i_actor, new%position, old)
    end subroutine Translation_update_actor

    subroutine Translation_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%translation%num_hits = changes_counters%translation%num_hits + 1
    end subroutine Translation_increment_hit

    subroutine Translation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%translation%num_successes = changes_counters%translation%num_successes + 1
    end subroutine Translation_increment_success

!end implementation Concrete_One_Particle_Translation

!implementation Concrete_One_Particle_Rotation

    subroutine Rotation_define_change(this, abort, new, old, i_actor)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_actor

        if (this%mixture%components(i_actor)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        old%i = random_integer(this%mixture%components(i_actor)%orientations%get_num())
        old%position = this%mixture%components(i_actor)%positions%get(old%i)
        old%orientation = this%mixture%components(i_actor)%orientations%get(old%i)
        old%dipole_moment = this%mixture%components(i_actor)%dipole_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%changes_components(i_actor)%rotated_orientations%get(new%i)
        new%dipole_moment = this%mixture%components(i_actor)%dipole_moments%get_norm() * new%&
            orientation
    end subroutine Rotation_define_change

    subroutine Rotation_visit_walls(this, overlap, delta_energy, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        delta_energy = 0._DP
    end subroutine Rotation_visit_walls

    subroutine Rotation_visit_short(this, overlap, delta_energies, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        delta_energies = 0._DP
    end subroutine Rotation_visit_short

    subroutine Rotation_visit_field(this, delta_energy, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta_energy = dipoles_field_visit_rotation(this%environment%external_field, new%&
            dipole_moment, old)
    end subroutine Rotation_visit_field

    subroutine Rotation_visit_dipolar(this, delta_energies, delta_shared_energy, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: real_energies_new, real_energies_old
        integer :: i_component, i_exclude

        do i_component = 1, size(this%dipolar_interactions_dynamic%real_components, 1)
            i_exclude = merge(new%i, 0, i_component == i_actor)
            call this%dipolar_interactions_dynamic%real_components(i_component, i_actor)%component%&
                visit(real_energies_new(i_component), new, visit_different, i_exclude)
            call this%dipolar_interactions_dynamic%real_components(i_component, i_actor)%component%&
                visit(real_energies_old(i_component), old, visit_different, i_exclude)
        end do
        delta_energies = real_energies_new - real_energies_old
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitor%&
                visit_rotation(i_actor, new%dipole_moment, old) + &
            this%dipolar_interactions_dynamic%surf_mixture%&
                visit_rotation(i_actor, new%dipole_moment, old%dipole_moment) - &
            this%dipolar_interactions_dynamic%dlc_visitor%&
                visit_rotation(i_actor, new%dipole_moment, old)
    end subroutine Rotation_visit_dipolar

    subroutine Rotation_update_actor(this, i_actor, new, old)
        class(Concrete_One_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_actor
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        call this%mixture%components(i_actor)%orientations%set(new%i, new%orientation)
        call this%mixture%total_moment%remove(i_actor, old%dipole_moment)
        call this%mixture%total_moment%add(i_actor, new%dipole_moment)
        call this%dipolar_interactions_static%reci_structure%&
            update_rotation(i_actor, new%dipole_moment, old)
        call this%dipolar_interactions_static%dlc_structures%&
            update_rotation(i_actor, new%dipole_moment, old)
    end subroutine Rotation_update_actor

    subroutine Rotation_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_hits = changes_counters%rotation%num_hits + 1
    end subroutine Rotation_increment_hit

    subroutine Rotation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_successes = changes_counters%rotation%num_successes + 1
    end subroutine Rotation_increment_success

!end implementation Concrete_One_Particle_Rotation

!implementation Null_One_Particle_Move

    subroutine Null_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes_components, can_move, &
        selector)
        class(Null_One_Particle_Move), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static
        type(Changes_Component_Wrapper), target, intent(in) :: changes_components(:)
        logical, intent(in) :: can_move(:)
        class(Abstract_Tower_Sampler), intent(in) :: selector
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_One_Particle_Move), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_reset_selector(this)
        class(Null_One_Particle_Move), intent(inout) :: this
    end subroutine Null_reset_selector

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_One_Particle_Move), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_One_Particle_Move), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

    subroutine Null_metropolis_algorithm(this, success, deltas, i_actor)
        class(Null_One_Particle_Move), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_actor
        success = .false.
        deltas%field_energy = 0._DP; deltas%walls_energy = 0._DP
        deltas%short_energies = 0._DP; deltas%dipolar_energies = 0._DP
        deltas%dipolar_shared_energy = 0._DP
    end subroutine Null_metropolis_algorithm

    subroutine Null_define_change(this, abort, new, old, i_actor)
         class(Null_One_Particle_Move), intent(in) :: this
         logical, intent(out) :: abort
         type(Concrete_Temporary_Particle), intent(out) :: new, old
         integer, intent(in) :: i_actor
         abort = .true.
         new%i = 0; new%position = 0._DP; new%orientation = 0._DP; new%dipole_moment = 0._DP
         old%i = 0; old%position = 0._DP; old%orientation = 0._DP; old%dipole_moment = 0._DP
     end subroutine Null_define_change

     subroutine Null_visit_walls(this, overlap, delta_energy, i_actor, new, old)
         class(Null_One_Particle_Move), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: delta_energy
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         overlap = .false.
         delta_energy = 0._DP
     end subroutine Null_visit_walls

     subroutine Null_visit_short(this, overlap, delta_energies, i_actor, new, old)
         class(Null_One_Particle_Move), intent(in) :: this
         logical, intent(out) :: overlap
         real(DP), intent(out) :: delta_energies(:)
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         overlap = .false.
         delta_energies = 0._DP
     end subroutine Null_visit_short

    subroutine Null_visit_field(this, delta_energy, new, old)
        class(Null_One_Particle_Move), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        delta_energy = 0._DP
    end subroutine Null_visit_field

     subroutine Null_visit_dipolar(this, delta_energies, delta_shared_energy, i_actor, new, old)
         class(Null_One_Particle_Move), intent(in) :: this
         real(DP), intent(out) :: delta_energies(:)
         real(DP), intent(out) :: delta_shared_energy
         type(Concrete_Temporary_Particle), intent(in) :: new, old
         integer, intent(in) :: i_actor
         delta_energies = 0._DP; delta_shared_energy = 0._DP
     end subroutine Null_visit_dipolar

     subroutine Null_update_actor(this, i_actor, new, old)
         class(Null_One_Particle_Move), intent(in) :: this
         integer, intent(in) :: i_actor
         type(Concrete_Temporary_Particle), intent(in) :: new, old
     end subroutine Null_update_actor

    subroutine Null_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_increment_hit

    subroutine Null_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters
    end subroutine Null_increment_success

!end implementation Null_One_Particle_Move

end module classes_one_particle_move
