module classes_box_particle_exchange

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_particle_wrapper, only: Concrete_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
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
use procedures_exchange_visitors, only: exchange_visit_add => visit_add, exchange_visit_remove => &
    visit_remove
use procedures_metropolis_algorithm, only: metropolis_algorithm
use procedures_exchange_updaters, only: exchange_update_add => update_add, &
    exchange_update_remove => update_remove

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Abstract_Box_Particle_Exchange
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic(:) => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        type(Changes_Wrapper), pointer :: changes => null()
        logical, allocatable :: can_exchange(:, :)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset_selectors => Abstract_reset_selectors
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure(Abstract_define_exchange), private, deferred :: define_exchange
        procedure(Abstract_acceptation_probability), private, deferred :: acceptation_probability
        procedure(Abstract_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_visit_short), private, deferred :: visit_short
        procedure(Abstract_visit_field), private, deferred :: visit_field
        procedure(Abstract_visit_dipolar), private, deferred :: visit_dipolar
        procedure(Abstract_update_component), private, deferred :: update_component
        procedure(Abstract_increment_hit), private, nopass, deferred :: increment_hit
        procedure(Abstract_increment_success), private, nopass, deferred :: increment_success
    end type Abstract_Box_Particle_Exchange

    abstract interface

        subroutine Abstract_define_exchange(this, abort, particle, i_box, i_component)
        import :: Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: abort
            type(Concrete_Particle), intent(out) :: particle
            integer, intent(in) :: i_box, i_component
        end subroutine Abstract_define_exchange

        pure real(DP) function Abstract_acceptation_probability(this, i_box, i_component, &
            delta_energy)
        import :: DP, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            integer, intent(in) :: i_box, i_component
            real(DP), intent(in) :: delta_energy
        end function Abstract_acceptation_probability

        pure subroutine Abstract_visit_walls(this, overlap, delta_energy, i_box, i_component, &
            particle)
        import :: DP, Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        import :: DP, Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_box, i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_visit_short

        pure subroutine Abstract_visit_field(this, delta_energy, i_box, particle)
        import :: DP, Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            real(DP), intent(out) :: delta_energy
            integer, intent(in) :: i_box
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_visit_field

        pure subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
            i_component, particle)
        import :: DP, Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:)
            real(DP), intent(out) :: delta_shared_energy
            integer, intent(in) :: i_box, i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_component(this, i_box, i_component, particle)
        import :: Concrete_Particle, Abstract_Box_Particle_Exchange
            class(Abstract_Box_Particle_Exchange), intent(in) :: this
            integer, intent(in) :: i_box, i_component
            type(Concrete_Particle), intent(in) :: particle
        end subroutine Abstract_update_component

        subroutine Abstract_increment_hit(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_increment_hit

        subroutine Abstract_increment_success(changes_counters)
        import :: Concrete_Changes_Counter
            type(Concrete_Changes_Counter), intent(inout) :: changes_counters
        end subroutine Abstract_increment_success

    end interface

    type, extends(Abstract_Box_Particle_Exchange), public :: Box_Particle_Add
    contains
        procedure, private :: define_exchange => Add_define_exchange
        procedure, private :: acceptation_probability => Add_acceptation_probability
        procedure, private :: visit_walls => Add_visit_walls
        procedure, private :: visit_short => Add_visit_short
        procedure, private :: visit_field => Add_visit_field
        procedure, private :: visit_dipolar => Add_visit_dipolar
        procedure, private :: update_component => Add_update_component
        procedure, nopass, private :: increment_hit => Add_increment_hit
        procedure, nopass, private :: increment_success => Add_increment_success
    end type Box_Particle_Add

    type, extends(Abstract_Box_Particle_Exchange), public :: Box_Particle_Remove
    contains
        procedure, private :: define_exchange => Remove_define_exchange
        procedure, private :: acceptation_probability => Remove_acceptation_probability
        procedure, private :: visit_walls => Remove_visit_walls
        procedure, private :: visit_short => Remove_visit_short
        procedure, private :: visit_field => Remove_visit_field
        procedure, private :: visit_dipolar => Remove_visit_dipolar
        procedure, private :: update_component => Remove_update_component
        procedure, nopass, private :: increment_hit => Remove_increment_hit
        procedure, nopass, private :: increment_success => Remove_increment_success
    end type Box_Particle_Remove

contains

!implementation Abstract_Box_Particle_Exchange

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes, can_exchange, selectors)
        class(Abstract_Box_Particle_Exchange), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic(:)
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
        class(Abstract_Box_Particle_Exchange), intent(inout) :: this

        call tower_sampler_destroy(this%selectors)
        if (allocated(this%can_exchange)) deallocate(this%can_exchange)
        this%changes => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selectors(this)
        class(Abstract_Box_Particle_Exchange), intent(inout) :: this

        call selectors_reset(this%selectors, this%mixture%components, this%can_exchange)
    end subroutine Abstract_reset_selectors

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Box_Particle_Exchange), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%selectors)
            num_choices = num_choices + this%selectors(i_box)%get_num_choices()
        end do
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Particle_Exchange), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Energies) :: deltas
        integer :: i_box, i_component

        i_box = random_integer(size(this%environment%periodic_boxes))
        i_component = this%selectors(i_box)%get()
        call this%increment_hit(observables%changes(i_box)%changes_counters(i_component))
        allocate(deltas%short_energies(size(observables%energies(i_box)%short_energies)))
        allocate(deltas%dipolar_energies(size(observables%energies(i_box)%dipolar_energies)))

        call this%metropolis_algorithm(success, deltas, i_box, i_component)

        if (success) then
            observables%nums_particles(i_component, i_box) = this%mixture%components(i_component, &
                i_box)%num_particles%get()
            call observables_energies_set(observables%energies(i_box), deltas, i_component)
            call this%increment_success(observables%changes(i_box)%changes_counters(i_component))
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, i_box, i_component)
        class(Abstract_Box_Particle_Exchange), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_box, i_component

        type(Concrete_Particle) :: particle
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
        success = metropolis_algorithm(this%acceptation_probability(i_box, i_component, &
            delta_energy))
        if (success) call this%update_component(i_box, i_component, particle)
    end subroutine Abstract_metropolis_algorithm

!end implementation Abstract_Box_Particle_Exchange

!implementation Box_Particle_Add

    subroutine Add_define_exchange(this, abort, particle, i_box, i_component)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Particle), intent(out) :: particle
        integer, intent(in) :: i_box, i_component

        abort = .false.
        particle%i = this%mixture%components(i_component, i_box)%num_particles%get() + 1
        particle%position = this%changes%random_positions(i_box)%get(i_component)
        particle%orientation = this%changes%random_orientation%get(i_component)
        particle%dipole_moment = this%mixture%components(i_component, i_box)%dipole_moments%&
            get_norm() * particle%orientation
    end subroutine Add_define_exchange

    !> \[
    !>      P[N \to N+1] = \min \left( 1,
    !>          \frac{V \rho}{N+1} e^{-\beta \Delta U_{N \to N+1}} a^{N} \right)
    !> \]
    pure real(DP) function Add_acceptation_probability(this, i_box, i_component, delta_energy) &
        result(probability)
        class(Box_Particle_Add), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        real(DP), intent(in) :: delta_energy

        associate(component => this%mixture%components(i_component, i_box))
            probability = product(this%environment%accessible_domains(i_box)%get_size()) * &
                component%chemical_potential%get_density() / &
                (real(component%num_particles%get() + 1, DP)) * &
                exp(-delta_energy / this%environment%temperature%get()) / &
                component%chemical_potential%get_inv_pow_activity()
        end associate
        probability = min(1._DP, probability)
    end function Add_acceptation_probability

    pure subroutine Add_visit_walls(this, overlap, delta_energy, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call this%environment%visitable_walls(i_box)%visit(overlap, delta_energy, particle%&
            position, this%short_interactions%wall_pairs(i_component)%potential)
    end subroutine Add_visit_walls

    subroutine Add_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_visit_add(overlap, delta_energies, i_component, particle, this%&
            short_interactions%cells(i_box))
    end subroutine Add_visit_short

    pure subroutine Add_visit_field(this, delta_energy, i_box, particle)
        class(Box_Particle_Add), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Particle), intent(in) :: particle

        delta_energy = dipoles_field_visit_add(this%environment%external_fields(i_box), particle)
    end subroutine Add_visit_field

    pure subroutine Add_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
        i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_visit_add(delta_energies, delta_shared_energy, i_component, particle, this%&
            dipolar_interactions_dynamic(i_box))
    end subroutine Add_visit_dipolar

    subroutine Add_update_component(this, i_box, i_component, particle)
        class(Box_Particle_Add), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_update_add(this%mixture%components(i_component, i_box), &
            this%mixture%total_moments(i_box), this%short_interactions%cells(i_box), &
            this%dipolar_interactions_static(i_box), i_component, particle)
    end subroutine Add_update_component

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
        type(Concrete_Particle), intent(out) :: particle
        logical, intent(out) :: abort
        integer, intent(in) :: i_box, i_component

        if (this%mixture%components(i_component, i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        particle%i = random_integer(this%mixture%components(i_component, i_box)%num_particles%get())
        particle%position = this%mixture%components(i_component, i_box)%positions%get(particle%i)
        particle%orientation = this%mixture%components(i_component, i_box)%orientations%&
            get(particle%i)
        particle%dipole_moment = this%mixture%components(i_component, i_box)%dipole_moments%&
            get(particle%i)
    end subroutine Remove_define_exchange

    !> \[
    !>      P[N \to N-1] = \min \left( 1,
    !>          \frac{N}{V \rho} e^{-\beta \Delta U_{N \to N-1}} a^{-N} \right)
    !> \]
    pure real(DP) function Remove_acceptation_probability(this, i_box, i_component, delta_energy) &
        result(probability)
        class(Box_Particle_Remove), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        real(DP), intent(in) :: delta_energy

        associate(component => this%mixture%components(i_component, i_box))
            probability = real(component%num_particles%get(), DP) / &
                product(this%environment%accessible_domains(i_box)%get_size()) / &
                component%chemical_potential%get_density() * &
                exp(-delta_energy / this%environment%temperature%get()) * &
                component%chemical_potential%get_inv_pow_activity()
        end associate
        probability = min(1._DP, probability)
    end function Remove_acceptation_probability

    pure subroutine Remove_visit_walls(this, overlap, delta_energy, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call this%environment%visitable_walls(i_box)%visit(overlap, delta_energy, particle%&
            position, this%short_interactions%wall_pairs(i_component)%potential)
        delta_energy = -delta_energy
    end subroutine Remove_visit_walls

    subroutine Remove_visit_short(this, overlap, delta_energies, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_visit_remove(overlap, delta_energies, i_component, particle, this%&
            short_interactions%cells(i_box))
    end subroutine Remove_visit_short

    pure subroutine Remove_visit_field(this, delta_energy, i_box, particle)
        class(Box_Particle_Remove), intent(in) :: this
        real(DP), intent(out) :: delta_energy
        integer, intent(in) :: i_box
        type(Concrete_Particle), intent(in) :: particle

        delta_energy = dipoles_field_visit_remove(this%environment%external_fields(i_box), particle)
    end subroutine Remove_visit_field

    pure subroutine Remove_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
        i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_visit_remove(delta_energies, delta_shared_energy, i_component, particle, &
            this%dipolar_interactions_dynamic(i_box))
    end subroutine Remove_visit_dipolar

    subroutine Remove_update_component(this, i_box, i_component, particle)
        class(Box_Particle_Remove), intent(in) :: this
        integer, intent(in) :: i_box, i_component
        type(Concrete_Particle), intent(in) :: particle

        call exchange_update_remove(this%mixture%components(i_component, i_box), &
            this%mixture%total_moments(i_box), this%short_interactions%cells(i_box), &
            this%dipolar_interactions_static(i_box), i_component, particle)
    end subroutine Remove_update_component

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
