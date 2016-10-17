module classes_box_particle_exchange

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
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add, &
    dipoles_field_visit_remove => visit_remove
use module_changes_success, only: Concrete_Changes_Counter
use types_changes_wrapper, only: Changes_Wrapper
use types_observables_energies, only: Concrete_Single_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Box_Particle_Exchange
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        type(Changes_Wrapper), pointer :: changes => null()
        logical, allocatable :: can_exchange(:, :)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure :: reset_selector => Abstract_reset_selector
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure(Abstract_define_exchange), private, deferred :: define_exchange
        procedure(Abstract_acceptation_probability), private, deferred :: acceptation_probability
        procedure(Abstract_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_visit_short), private, deferred :: visit_short
        procedure(Abstract_visit_field), private, deferred :: visit_field
        procedure(Abstract_visit_dipolar), private, deferred :: visit_dipolar
        procedure(Abstract_update_actor), private, deferred :: update_actor
        procedure(Abstract_increment_hit), private, nopass, deferred :: increment_hit
        procedure(Abstract_increment_success), private, nopass, deferred :: increment_success
    end type Box_Particle_Exchange

    abstract interface

        subroutine Abstract_define_exchange(this, abort, particle, i_box, i_component)
        import :: Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: abort
            type(Concrete_Temporary_Particle), intent(out) :: particle
            integer, intent(in) :: i_box, i_component
        end subroutine Abstract_define_exchange

        pure real(DP) function Abstract_acceptation_probability(this, i_box, i_component, delta_energy)
        import :: DP, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            integer, intent(in) :: i_box, i_component
            real(DP), intent(in) :: delta_energy
        end function Abstract_acceptation_probability

        subroutine Abstract_visit_walls(this, overlap, delta_energy, i_box, i_component, particle)
        import :: DP, Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: particle
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        import :: DP, Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: particle
        end subroutine Abstract_visit_short

        subroutine Abstract_visit_field(this, delta_energy, i_box, particle)
        import :: DP, Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box
            type(Concrete_Temporary_Particle), intent(in) :: particle
        end subroutine Abstract_visit_field

        subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, i_component, &
            particle)
        import :: DP, Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:)
            real(DP), intent(out) :: delta_shared_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: particle
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_actor(this, i_box, i_component, particle)
        import :: Concrete_Temporary_Particle, Box_Particle_Exchange
            class(Box_Particle_Exchange), intent(in) :: this
            integer, intent(in) :: i_box, i_component
            type(Concrete_Temporary_Particle), intent(in) :: particle
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

    type, extends(Box_Particle_Exchange), public :: Box_Particle_Add
    contains
        procedure, private :: define_exchange => Add_define_exchange
        procedure, private :: acceptation_probability => Add_acceptation_probability
        procedure, private :: visit_walls => Add_visit_walls
        procedure, private :: visit_short => Add_visit_short
        procedure, private :: visit_field => Add_visit_field
        procedure, private :: visit_dipolar => Add_visit_dipolar
        procedure, private :: update_actor => Add_update_actor
        procedure, nopass, private :: increment_hit => Add_increment_hit
        procedure, nopass, private :: increment_success => Add_increment_success
    end type Box_Particle_Add

    type, extends(Box_Particle_Exchange), public :: Box_Particle_Remove
    contains
        procedure, private :: define_exchange => Remove_define_exchange
        procedure, private :: acceptation_probability => Remove_acceptation_probability
        procedure, private :: visit_walls => Remove_visit_walls
        procedure, private :: visit_short => Remove_visit_short
        procedure, private :: visit_field => Remove_visit_field
        procedure, private :: visit_dipolar => Remove_visit_dipolar
        procedure, private :: update_actor => Remove_update_actor
        procedure, nopass, private :: increment_hit => Remove_increment_hit
        procedure, nopass, private :: increment_success => Remove_increment_success
    end type Box_Particle_Remove

contains

!implementation Box_Particle_Exchange

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes, can_exchange, &
        selectors)
        class(Box_Particle_Exchange), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: &
            dipolar_interactions_static(:)
        type(Changes_Wrapper), target, intent(in) :: changes
        logical, intent(in) :: can_exchange(:, :)
        class(Abstract_Tower_Sampler), intent(in) :: selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%changes => changes
        allocate(this%can_exchange, source=can_exchange)
        allocate(this%selectors, source=selectors)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Box_Particle_Exchange), intent(inout) :: this

        call tower_sampler_destroy(this%selectors)
        if (allocated(this%can_exchange)) deallocate(this%can_exchange)
        this%changes => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selector(this)
        class(Box_Particle_Exchange), intent(inout) :: this

        call selectors_reset(this%selectors, this%mixture%gemc_components, this%can_exchange)
    end subroutine Abstract_reset_selector

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Box_Particle_Exchange), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%selectors)
            num_choices = num_choices + this%selectors(i_box)%get_num_choices()
        end do
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Box_Particle_Exchange), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Energies) :: deltas
        integer :: i_box, i_component

        i_box = random_integer(size(this%environment%periodic_boxes))
        i_component = this%selectors(i_box)%get()
        call this%increment_hit(observables%changes_counters(i_component))
        allocate(deltas%short_energies(size(observables%energies%short_energies)))
        allocate(deltas%dipolar_energies(size(observables%energies%dipolar_energies)))
        call this%metropolis_algorithm(success, deltas, i_box, i_component)
        if (success) then
            observables%nums_particles(i_component) = this%mixture%gemc_components(i_component, i_box)%&
                num_particles%get()
            call observables_energies_set(observables%energies, deltas, i_component)
            call this%increment_success(observables%changes_counters(i_component))
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, i_box, i_component)
        class(Box_Particle_Exchange), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_box, i_component

        type(Concrete_Temporary_Particle) :: particle
        real(DP) :: delta_energy
        logical :: abort, overlap

        success = .false.
        call this%define_exchange(abort, particle, i_box, i_component)
        if (abort) return

        call this%visit_walls(overlap, deltas%walls_energy, i_box, i_component, particle)
        if (overlap) return
        call this%visit_short(overlap, deltas%short_energies, i_box, i_component, particle)
        if (overlap) return
        call this%visit_field(deltas%field_energy, i_box, particle)
        call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, i_box, &
            i_component, particle)

        delta_energy = deltas%walls_energy + deltas%field_energy + &
            sum(deltas%short_energies + deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(this%acceptation_probability(i_box, i_component, delta_energy))
        if (success) call this%update_actor(i_box, i_component, particle)
    end subroutine Abstract_metropolis_algorithm

!end implementation Box_Particle_Exchange

!implementation Box_Particle_Add

    subroutine Add_define_exchange(this, abort, particle, i_box, i_component)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: particle
        integer, intent(in) :: i_box, i_component

        abort = .false.
        particle%i = this%mixture%gemc_components(i_component, i_box)%num_particles%get() + 1
        particle%position = this%changes%random_positions(i_box)%get(i_component)
        particle%orientation = this%changes%random_orientation%get(i_component)
        particle%dipole_moment = this%mixture%gemc_components(i_component, i_box)%dipole_moments%get_norm()*&
            particle%orientation
    end subroutine Add_define_exchange

    !> \[
    !>      P[N \to N+1] = min \left( 1,
    !>          \frac{V \rho}{N+1} e^{-\beta \Delta U_{N \to N+1}} a^{N} \right)
    !> \]
    pure real(DP) function Add_acceptation_probability(this, i_box, i_component, delta_energy) &
        result(probability)
        class(Box_Particle_Add), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        real(DP), intent(in) :: delta_energy

        associate(component => this%mixture%gemc_components(i_component, i_box))
            probability = product(this%environment%accessible_domains(i_box)%get_size()) * &
                component%chemical_potential%get_density() / &
                (real(component%num_particles%get() + 1, DP)) * &
                exp(-delta_energy / this%environment%temperature%get()) / &
                component%chemical_potential%get_inv_pow_activity()
        end associate
        probability = min(1._DP, probability)
    end function Add_acceptation_probability

    subroutine Add_visit_walls(this, overlap, delta_energy, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%environment%gemc_visitable_walls(i_box)%visit(overlap, delta_energy, particle%&
            position, this%short_interactions%wall_pairs(i_component)%potential)
    end subroutine Add_visit_walls

    subroutine Add_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component, i_exclude

        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call this%short_interactions%cells(i_box)%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, delta_energies(j_component), particle, visit_different, &
                    i_exclude)
            if (overlap) return
        end do
    end subroutine Add_visit_short

    subroutine Add_visit_field(this, delta_energy, i_box, particle)
        class(Box_Particle_Add), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: particle

        delta_energy = dipoles_field_visit_add(this%environment%external_fields(i_box), particle)
    end subroutine Add_visit_field

    subroutine Add_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions_dynamic%gemc_real_components, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call this%dipolar_interactions_dynamic%gemc_real_components(j_component, i_component, i_box)%component%&
                visit(delta_energies(j_component), particle, visit_different, i_exclude)
        end do
        delta_energies(i_component) = delta_energies(i_component) - this%dipolar_interactions_dynamic%&
            gemc_self_components(i_component, i_box)%component%meet(particle%dipole_moment)
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitors(i_box)%visit_add(i_component, particle) + &
            this%dipolar_interactions_dynamic%gemc_surf_mixture(i_box)%visit_add(i_component, particle%&
                dipole_moment) - &
            this%dipolar_interactions_dynamic%dlc_visitors(i_box)%visit_add(i_component, particle)
    end subroutine Add_visit_dipolar

    subroutine Add_update_actor(this, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component

        call this%mixture%gemc_components(i_component, i_box)%num_particles%set(this%mixture%&
            gemc_components(i_component, i_box)%num_particles%get() + 1)
        call this%mixture%gemc_components(i_component, i_box)%positions%add(particle%position)
        call this%mixture%gemc_components(i_component, i_box)%orientations%add(particle%orientation)
        call this%mixture%total_moments(i_box)%add(i_component, particle%dipole_moment)
        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 2)
            call this%short_interactions%cells(i_box)%visitable_cells(i_component, j_component)%add(particle)
        end do
        call this%dipolar_interactions_static(i_box)%reci_structure%update_add(i_component, particle)
        call this%dipolar_interactions_static(i_box)%dlc_structures%update_add(i_component, particle)
    end subroutine Add_update_actor

    subroutine Add_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%add%num_hits = changes_counters%add%num_hits + 1
    end subroutine Add_increment_hit

    subroutine Add_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%add%num_successes = changes_counters%add%num_successes + 1
    end subroutine Add_increment_success

!end implementation Box_Particle_Add

!implementation Box_Particle_Remove

    subroutine Remove_define_exchange(this, abort, particle, i_box, i_component)
        class(Box_Particle_Remove), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(out) :: particle
        logical, intent(out) :: abort
        integer, intent(in) :: i_box, i_component

        if (this%mixture%gemc_components(i_component, i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        particle%i = random_integer(this%mixture%gemc_components(i_component, i_box)%num_particles%get())
        particle%position = this%mixture%gemc_components(i_component, i_box)%positions%get(particle%i)
        particle%orientation = this%mixture%gemc_components(i_component, i_box)%orientations%get(particle%i)
        particle%dipole_moment = this%mixture%gemc_components(i_component, i_box)%dipole_moments%&
            get(particle%i)
    end subroutine Remove_define_exchange

    !> \[
    !>      P[N \to N-1] = min \left( 1,
    !>          \frac{N}{V \rho} e^{-\beta \Delta U_{N \to N-1}} a^{-N} \right)
    !> \]
    pure real(DP) function Remove_acceptation_probability(this, i_box, i_component, delta_energy) &
        result(probability)
        class(Box_Particle_Remove), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        real(DP), intent(in) :: delta_energy

        associate(component => this%mixture%gemc_components(i_component, i_box))
            probability = real(component%num_particles%get(), DP) / &
                product(this%environment%accessible_domains(i_box)%get_size()) / &
                component%chemical_potential%get_density() * &
                exp(-delta_energy / this%environment%temperature%get()) * &
                component%chemical_potential%get_inv_pow_activity()
        end associate
        probability = min(1._DP, probability)
    end function Remove_acceptation_probability

    subroutine Remove_visit_walls(this, overlap, delta_energy, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%environment%gemc_visitable_walls(i_box)%visit(overlap, delta_energy, particle%&
            position, this%short_interactions%wall_pairs(i_component)%potential)
        delta_energy = -delta_energy
    end subroutine Remove_visit_walls

    subroutine Remove_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component, i_exclude

        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call this%short_interactions%cells(i_box)%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, delta_energies(j_component), particle, visit_different, &
                    i_exclude)
        end do
        delta_energies = -delta_energies
    end subroutine Remove_visit_short

    subroutine Remove_visit_field(this, delta_energy, i_box, particle)
        class(Box_Particle_Remove), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: particle

        delta_energy = dipoles_field_visit_remove(this%environment%external_fields(i_box), particle)
    end subroutine Remove_visit_field

    subroutine Remove_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions_dynamic%gemc_real_components, 1)
            i_exclude = merge(particle%i, 0, j_component == i_component)
            call this%dipolar_interactions_dynamic%gemc_real_components(j_component, i_component, i_box)%component%&
                visit(delta_energies(j_component), particle, visit_different, i_exclude)
        end do
        delta_energies(i_component) = delta_energies(i_component) - this%dipolar_interactions_dynamic%&
            gemc_self_components(i_component, i_box)%component%meet(particle%dipole_moment)
        delta_energies = -delta_energies
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitors(i_box)%visit_remove(i_component, particle) + &
            this%dipolar_interactions_dynamic%gemc_surf_mixture(i_box)%visit_remove(i_component, particle%&
                dipole_moment) - &
            this%dipolar_interactions_dynamic%dlc_visitors(i_box)%visit_remove(i_component, particle)
    end subroutine Remove_visit_dipolar

    subroutine Remove_update_actor(this, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        integer :: j_component

        call this%dipolar_interactions_static(i_box)%dlc_structures%update_remove(i_component, particle)
        call this%dipolar_interactions_static(i_box)%reci_structure%update_remove(i_component, particle)
        do j_component = size(this%short_interactions%cells(i_box)%visitable_cells, 2), 1, -1
            call this%short_interactions%cells(i_box)%visitable_cells(i_component, j_component)%remove(particle)
        end do
        call this%mixture%total_moments(i_box)%remove(i_component, particle%dipole_moment)
        call this%mixture%gemc_components(i_component, i_box)%orientations%remove(particle%i)
        call this%mixture%gemc_components(i_component, i_box)%positions%remove(particle%i)
        call this%mixture%gemc_components(i_component, i_box)%num_particles%set(this%mixture%&
            gemc_components(i_component, i_box)%num_particles%get() - 1)
    end subroutine Remove_update_actor

    subroutine Remove_increment_hit(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%remove%num_hits = changes_counters%remove%num_hits + 1
    end subroutine Remove_increment_hit

    subroutine Remove_increment_success(changes_counters)
        type(Concrete_Changes_Counter), intent(inout) :: changes_counters

        changes_counters%remove%num_successes = changes_counters%remove%num_successes + 1
    end subroutine Remove_increment_success

!end implementation Box_Particle_Remove

end module classes_box_particle_exchange
