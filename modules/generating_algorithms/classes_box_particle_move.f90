module classes_box_particle_move

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
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Abstract_Box_Particle_Move
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic(:) => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        type(Changes_Component_Wrapper), pointer :: changes_components(:, :) => null()
        logical, allocatable :: can_move(:, :)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure :: reset_selectors => Abstract_reset_selectors
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
    end type Abstract_Box_Particle_Move

    abstract interface

        subroutine Abstract_define_change(this, abort, new, old, i_box, i_component)
        import :: Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            logical, intent(out) :: abort
            type(Concrete_Temporary_Particle), intent(out) :: new, old
            integer, intent(in) :: i_box, i_component
        end subroutine Abstract_define_change

        subroutine Abstract_visit_walls(this, overlap, delta_energy, i_box, i_component, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, delta_energies, i_box, i_component, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_short

        subroutine Abstract_visit_field(this, delta_energy, i_box, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_field

        subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
            i_component, new, old)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:)
            real(DP), intent(out) :: delta_shared_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: new, old
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_actor(this, i_box, i_component, new, old)
        import :: Concrete_Temporary_Particle, Abstract_Box_Particle_Move
            class(Abstract_Box_Particle_Move), intent(in) :: this
            integer, intent(in) :: i_box, i_component
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

    type, extends(Abstract_Box_Particle_Move), public :: Box_Particle_Translation
    contains
        procedure, private :: define_change => Translation_define_change
        procedure, private :: visit_walls => Translation_visit_walls
        procedure, private :: visit_short => Translation_visit_short
        procedure, private :: visit_field => Translation_visit_field
        procedure, private :: visit_dipolar => Translation_visit_dipolar
        procedure, private :: update_actor => Translation_update_actor
        procedure, private, nopass :: increment_hit => Translation_increment_hit
        procedure, private, nopass :: increment_success => Translation_increment_success
    end type Box_Particle_Translation

    type, extends(Abstract_Box_Particle_Move), public :: Box_Particle_Rotation
    contains
        procedure, private :: define_change => Rotation_define_change
        procedure, private :: visit_walls => Rotation_visit_walls
        procedure, private :: visit_short => Rotation_visit_short
        procedure, private :: visit_field => Rotation_visit_field
        procedure, private :: visit_dipolar => Rotation_visit_dipolar
        procedure, private :: update_actor => Rotation_update_actor
        procedure, private, nopass :: increment_hit => Rotation_increment_hit
        procedure, private, nopass :: increment_success => Rotation_increment_success
    end type Box_Particle_Rotation

contains

!implementation Abstract_Box_Particle_Move

    !> @note this%selectors construction is delayed in
    !> [[classes_box_particle_move:Abstract_reset_selectors]]
    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes_components, can_move, &
        selectors)
        class(Abstract_Box_Particle_Move), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic(:)
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: &
            dipolar_interactions_static(:)
        type(Changes_Component_Wrapper), target, intent(in) :: changes_components(:, :)
        logical, intent(in) :: can_move(:, :)
        class(Abstract_Tower_Sampler), intent(in) :: selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%changes_components => changes_components
        allocate(this%can_move, source=can_move)
        allocate(this%selectors, source=selectors)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Particle_Move), intent(inout) :: this

        call tower_sampler_destroy(this%selectors)
        if (allocated(this%can_move)) deallocate(this%can_move)
        this%changes_components => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selectors(this)
        class(Abstract_Box_Particle_Move), intent(inout) :: this

        call selectors_reset(this%selectors, this%mixture%gemc_components, this%can_move)
    end subroutine Abstract_reset_selectors

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Box_Particle_Move), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%selectors)
            num_choices = num_choices + this%selectors(i_box)%get_num_choices()
        end do
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Particle_Move), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Energies) :: deltas
        integer :: i_box, i_component

        i_box = random_integer(size(this%environment%periodic_boxes))
        i_component = this%selectors(i_box)%get()
        call this%increment_hit(observables%changes(i_box)%changes_counters(i_component))
        allocate(deltas%short_energies(size(observables%gemc_energies(i_box)%short_energies)))
        allocate(deltas%dipolar_energies(size(observables%gemc_energies(i_box)%dipolar_energies)))
        call this%metropolis_algorithm(success, deltas, i_box, i_component)
        if (success) then
            call observables_energies_set(observables%gemc_energies(i_box), deltas, i_component)
            call this%increment_success(observables%changes(i_box)%changes_counters(i_component))
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, i_box, i_component)
        class(Abstract_Box_Particle_Move), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_box, i_component

        real(DP) :: delta_energy
        type(Concrete_Temporary_Particle) :: new, old
        logical :: abort, overlap

        success = .false.
        call this%define_change(abort, new, old, i_box, i_component)
        if (abort) return

        call this%visit_walls(overlap, deltas%walls_energy, i_box, i_component, new, old)
        if (overlap) return
        call this%visit_short(overlap, deltas%short_energies, i_box, i_component, new, old)
        if (overlap) return
        call this%visit_field(deltas%field_energy, i_box, new, old)
        call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, i_box, &
            i_component, new, old)

        delta_energy = deltas%walls_energy + deltas%field_energy + &
            sum(deltas%short_energies + deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(min(1._DP, &
            exp(-delta_energy/this%environment%temperature%get())))
        if (success) call this%update_actor(i_box, i_component, new, old)
    end subroutine Abstract_metropolis_algorithm

!end implementation Abstract_Box_Particle_Move

!implementation Box_Particle_Translation

    subroutine Translation_define_change(this, abort, new, old, i_box, i_component)
        class(Box_Particle_Translation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_box, i_component

        if (this%mixture%gemc_components(i_component, i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        old%i = random_integer(this%mixture%gemc_components(i_component, i_box)%num_particles%get())
        old%position = this%mixture%gemc_components(i_component, i_box)%positions%get(old%i)
        old%orientation = this%mixture%gemc_components(i_component, i_box)%orientations%get(old%i)
        old%dipole_moment = this%mixture%gemc_components(i_component, i_box)%dipole_moments%get(old%i)
        new%i = old%i
        new%position = this%changes_components(i_component, i_box)%translated_positions%get(new%i)
        new%orientation = old%orientation
        new%dipole_moment = old%dipole_moment
    end subroutine Translation_define_change

    subroutine Translation_visit_walls(this, overlap, delta_energy, i_box, i_component, new, old)
        class(Box_Particle_Translation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: energy_new, energy_old

        call this%environment%gemc_visitable_walls(i_box)%visit(overlap, energy_new, new%position, this%&
            short_interactions%wall_pairs(i_component)%potential)
        if (overlap) return
        call this%environment%gemc_visitable_walls(i_box)%visit(overlap, energy_old, old%position, this%&
            short_interactions%wall_pairs(i_component)%potential)
        delta_energy = energy_new - energy_old
    end subroutine Translation_visit_walls

    subroutine Translation_visit_short(this, overlap, delta_energies, i_box, i_component, new, old)
        class(Box_Particle_Translation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: energies_new, energies_old
        integer :: j_component, i_exclude

        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
            i_exclude = merge(new%i, 0, j_component == i_component)
            call this%short_interactions%cells(i_box)%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, energies_new(j_component), new, visit_different, i_exclude)
            if (overlap) return
        end do
        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
            i_exclude = merge(old%i, 0, j_component == i_component)
            call this%short_interactions%cells(i_box)%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, energies_old(j_component), old, visit_different, i_exclude)
        end do
        delta_energies = energies_new - energies_old
    end subroutine Translation_visit_short

    subroutine Translation_visit_field(this, delta_energy, i_box, new, old)
        class(Box_Particle_Translation), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta_energy = dipoles_field_visit_translation(this%environment%external_fields(i_box), &
            new%position, old)
    end subroutine Translation_visit_field

    subroutine Translation_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, i_component,&
        new, old)
        class(Box_Particle_Translation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: real_energies_new, real_energies_old
        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions_dynamic(i_box)%real_components, 1)
            i_exclude = merge(new%i, 0, j_component == i_component)
            call this%dipolar_interactions_dynamic(i_box)%real_components(j_component, i_component)%&
                component%visit(real_energies_new(j_component), new, visit_different, i_exclude)
            call this%dipolar_interactions_dynamic(i_box)%real_components(j_component, i_component)%&
                component%visit(real_energies_old(j_component), old, visit_different, i_exclude)
        end do
        delta_shared_energy = &
            this%dipolar_interactions_dynamic(i_box)%reci_visitor%&
                visit_translation(i_component, new%position, old) - &
            this%dipolar_interactions_dynamic(i_box)%dlc_visitor%&
                visit_translation(i_component, new%position, old)
        delta_energies = real_energies_new - real_energies_old
    end subroutine Translation_visit_dipolar

    subroutine Translation_update_actor(this, i_box, i_component, new, old)
        class(Box_Particle_Translation), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        integer :: j_component

        call this%mixture%gemc_components(i_component, i_box)%positions%set(new%i, new%position)
        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 2)
            call this%short_interactions%cells(i_box)%visitable_cells(i_component, j_component)%&
                translate(new%position, old)
        end do
        call this%dipolar_interactions_static(i_box)%reci_structure%&
            update_translation(i_component, new%position, old)
        call this%dipolar_interactions_static(i_box)%dlc_structures%&
            update_translation(i_component, new%position, old)
    end subroutine Translation_update_actor

    subroutine Translation_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%translation%num_hits = changes_counters%translation%num_hits + 1
    end subroutine Translation_increment_hit

    subroutine Translation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%translation%num_successes = changes_counters%translation%num_successes + 1
    end subroutine Translation_increment_success

!end implementation Box_Particle_Translation

!implementation Box_Particle_Rotation

    subroutine Rotation_define_change(this, abort, new, old, i_box, i_component)
        class(Box_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new, old
        integer, intent(in) :: i_box, i_component

        if (this%mixture%gemc_components(i_component, i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        old%i = random_integer(this%mixture%gemc_components(i_component, i_box)%orientations%get_num())
        old%position = this%mixture%gemc_components(i_component, i_box)%positions%get(old%i)
        old%orientation = this%mixture%gemc_components(i_component, i_box)%orientations%get(old%i)
        old%dipole_moment = this%mixture%gemc_components(i_component, i_box)%dipole_moments%get(old%i)
        new%i = old%i
        new%position = old%position
        new%orientation = this%changes_components(i_component, i_box)%rotated_orientations%get(new%i)
        new%dipole_moment = this%mixture%gemc_components(i_component, i_box)%dipole_moments%get_norm() * &
            new%orientation
    end subroutine Rotation_define_change

    subroutine Rotation_visit_walls(this, overlap, delta_energy, i_box, i_component, new, old)
        class(Box_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        delta_energy = 0._DP
    end subroutine Rotation_visit_walls

    subroutine Rotation_visit_short(this, overlap, delta_energies, i_box, i_component, new, old)
        class(Box_Particle_Rotation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        overlap = .false.
        delta_energies = 0._DP
    end subroutine Rotation_visit_short

    subroutine Rotation_visit_field(this, delta_energy, i_box, new, old)
        class(Box_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        delta_energy = dipoles_field_visit_rotation(this%environment%external_fields(i_box), new%&
            dipole_moment, old)
    end subroutine Rotation_visit_field

    subroutine Rotation_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, i_component, &
        new, old)
        class(Box_Particle_Rotation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP), dimension(size(delta_energies)) :: real_energies_new, real_energies_old
        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions_dynamic(i_box)%real_components, 1)
            i_exclude = merge(new%i, 0, j_component == i_component)
            call this%dipolar_interactions_dynamic(i_box)%real_components(j_component, i_component)%&
                component%visit(real_energies_new(j_component), new, visit_different, i_exclude)
            call this%dipolar_interactions_dynamic(i_box)%real_components(j_component, i_component)%&
                component%visit(real_energies_old(j_component), old, visit_different, i_exclude)
        end do
        delta_energies = real_energies_new - real_energies_old
        delta_shared_energy = &
            this%dipolar_interactions_dynamic(i_box)%reci_visitor%&
                visit_rotation(i_component, new%dipole_moment, old) + &
            this%dipolar_interactions_dynamic(i_box)%surf_mixture%&
                visit_rotation(i_component, new%dipole_moment, old%dipole_moment) - &
            this%dipolar_interactions_dynamic(i_box)%dlc_visitor%&
                visit_rotation(i_component, new%dipole_moment, old)
    end subroutine Rotation_visit_dipolar

    subroutine Rotation_update_actor(this, i_box, i_component, new, old)
        class(Box_Particle_Rotation), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        call this%mixture%gemc_components(i_component, i_box)%orientations%set(new%i, new%orientation)
        call this%mixture%total_moments(i_box)%remove(i_component, old%dipole_moment)
        call this%mixture%total_moments(i_box)%add(i_component, new%dipole_moment)
        call this%dipolar_interactions_static(i_box)%reci_structure%&
            update_rotation(i_component, new%dipole_moment, old)
        call this%dipolar_interactions_static(i_box)%dlc_structures%&
            update_rotation(i_component, new%dipole_moment, old)
    end subroutine Rotation_update_actor

    subroutine Rotation_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_hits = changes_counters%rotation%num_hits + 1
    end subroutine Rotation_increment_hit

    subroutine Rotation_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%rotation%num_successes = changes_counters%rotation%num_successes + 1
    end subroutine Rotation_increment_success

!end implementation Box_Particle_Rotation

end module classes_box_particle_move
