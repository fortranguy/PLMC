module classes_box_particles_swap

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_destroy => destroy
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
use procedures_dipoles_field_interaction, only: dipoles_field_visit_translation => &
    visit_translation, dipoles_field_visit_add => visit_add, dipoles_field_visit_remove => &
    visit_remove
use types_changes_wrapper, only: Changes_Wrapper
use types_observables_energies, only: Concrete_Double_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use types_observables_changes, only: Concrete_Observables_Changes
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: Abstract_Box_Particles_Swap
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic(:) => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        type(Changes_Wrapper), pointer :: changes => null()
        logical, allocatable :: can_swap(:, :)
        class(Abstract_Hetero_Couples), allocatable :: couples(:)
        class(Abstract_Tower_Sampler), allocatable :: selectors(:) !! [i, j] <-> k: convert
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset_selectors => Abstract_reset_selectors
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure(Abstract_define_swap), private, deferred :: define_swap
        procedure(Abstract_acceptation_probability), private, deferred :: acceptation_probability
        procedure(Abstract_visit_walls), private, deferred :: visit_walls
        procedure(Abstract_visit_short), private, deferred :: visit_short
        procedure(Abstract_visit_field), private, deferred :: visit_field
        procedure(Abstract_visit_dipolar), private, deferred :: visit_dipolar
        procedure(Abstract_update_components), private, deferred :: update_components
        procedure(Abstract_increment_hit), private, deferred, nopass :: increment_hit
        procedure(Abstract_increment_success), private, deferred, nopass :: increment_success
    end type Abstract_Box_Particles_Swap

    abstract interface

        subroutine Abstract_define_swap(this, abort, particles, i_box, ij_components)
        import :: Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            logical, intent(out) :: abort
            type(Concrete_Temporary_Particle), intent(out) :: particles(:)
            integer, intent(in) :: i_box, ij_components(:)
        end subroutine Abstract_define_swap

        pure real(DP) function Abstract_acceptation_probability(this, i_box, ij_components, &
            delta_energy) result(probability)
        import :: DP, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            integer, intent(in) :: i_box, ij_components(:)
            real(DP), intent(in) :: delta_energy
        end function Abstract_acceptation_probability

        subroutine Abstract_visit_walls(this, overlap, delta_energies, i_box, ij_components, &
            particles)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_box, ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end subroutine Abstract_visit_walls

        subroutine Abstract_visit_short(this, overlap, delta_energies, i_box, ij_components, &
            particles)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: delta_energies(:, :)
            integer, intent(in) :: i_box, ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end subroutine Abstract_visit_short

        subroutine Abstract_visit_field(this, delta_energies, i_box, particles)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:)
            integer, intent(in) :: i_box
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end subroutine Abstract_visit_field

        subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
            ij_components, particles)
        import :: DP, Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            real(DP), intent(out) :: delta_energies(:, :)
            real(DP), intent(out) :: delta_shared_energy
            integer, intent(in) :: i_box, ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        end subroutine Abstract_visit_dipolar

        subroutine Abstract_update_components(this, i_box, ij_components, particles)
        import :: Concrete_Temporary_Particle, Abstract_Box_Particles_Swap
            class(Abstract_Box_Particles_Swap), intent(in) :: this
            integer, intent(in) :: i_box, ij_components(:)
            type(Concrete_Temporary_Particle), intent(in) ::  particles(:)
        end subroutine Abstract_update_components

        subroutine Abstract_increment_hit(changes, ij_components)
        import :: Concrete_Observables_Changes
            type(Concrete_Observables_Changes), intent(inout) :: changes
            integer, intent(in) :: ij_components(:)
        end subroutine Abstract_increment_hit

        subroutine Abstract_increment_success(changes, ij_components)
        import :: Concrete_Observables_Changes
            type(Concrete_Observables_Changes), intent(inout) :: changes
            integer, intent(in) :: ij_components(:)
        end subroutine Abstract_increment_success

    end interface

    !> Identities swap
    type, extends(Abstract_Box_Particles_Swap), public :: Box_Particles_Transmutation
    contains
        procedure, private :: define_swap => Transmutation_define_swap
        procedure, private :: acceptation_probability => Transmutation_acceptation_probability
        procedure, private :: visit_walls => Transmutation_visit_walls
        procedure, private :: visit_short => Transmutation_visit_short
        procedure, private :: visit_field => Transmutation_visit_field
        procedure, private :: visit_dipolar => Transmutation_visit_dipolar
        procedure, private :: update_components => Transmutation_update_components
        procedure, private, nopass :: increment_hit => Transmutation_increment_hit
        procedure, private, nopass :: increment_success => Transmutation_increment_success
    end type Box_Particles_Transmutation

    !> Positions swap
    type, extends(Abstract_Box_Particles_Swap), public :: Box_Particles_Switch
    contains
        procedure, private :: define_swap => Switch_define_swap
        procedure, private :: acceptation_probability => Switch_acceptation_probability
        procedure, private :: visit_walls => Switch_visit_walls
        procedure, private :: visit_short => Switch_visit_short
        procedure, private :: visit_field => Switch_visit_field
        procedure, private :: visit_dipolar => Switch_visit_dipolar
        procedure, private :: update_components => Switch_update_components
        procedure, private, nopass :: increment_hit => Switch_increment_hit
        procedure, private, nopass :: increment_success => Switch_increment_success
    end type Box_Particles_Switch

contains

!implementation Abstract_Box_Particles_Swap

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes, can_swap, couples, &
        selectors)
        class(Abstract_Box_Particles_Swap), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic(:)
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: &
            dipolar_interactions_static(:)
        type(Changes_Wrapper), target, intent(in) :: changes
        logical, intent(in) :: can_swap(:, :)
        class(Abstract_Hetero_Couples), intent(in) :: couples(:)
        class(Abstract_Tower_Sampler), intent(in) :: selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%changes => changes
        allocate(this%can_swap, source=can_swap)
        allocate(this%couples, source=couples)
        allocate(this%selectors, source=selectors)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Particles_Swap), intent(inout) :: this

        call tower_sampler_destroy(this%selectors)
        call hetero_couples_destroy(this%couples)
        if (allocated(this%can_swap)) deallocate(this%can_swap)
        this%changes => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_reset_selectors(this)
        class(Abstract_Box_Particles_Swap), intent(inout) :: this

        call selectors_reset(this%selectors, this%couples, this%mixture%components, this%can_swap)
    end subroutine Abstract_reset_selectors

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Box_Particles_Swap), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%selectors)
            num_choices = num_choices + this%selectors(i_box)%get_num_choices()
        end do
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Particles_Swap), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Double_Energies) :: deltas
        integer :: i_box, ij_components(2)

        i_box = random_integer(size(this%environment%periodic_boxes))
        ij_components = this%couples(i_box)%get(this%selectors(i_box)%get())
        call this%increment_hit(observables%changes(i_box), ij_components)
        allocate(deltas%short_energies(size(observables%energies(i_box)%short_energies), 2))
        allocate(deltas%dipolar_energies(size(observables%energies(i_box)%dipolar_energies), 2))
        call this%metropolis_algorithm(success, deltas, i_box, ij_components)
        if (success) then
            call observables_energies_set(observables%energies(i_box), deltas, ij_components)
            call this%increment_success(observables%changes(i_box), ij_components)
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, i_box, ij_components)
        class(Abstract_Box_Particles_Swap), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Energies), intent(inout) :: deltas
        integer, intent(in) :: i_box, ij_components(:)

        real(DP) :: delta_energy
        type(Concrete_Temporary_Particle) :: particles(2)
        logical :: abort, overlap

        success = .false.
        call this%define_swap(abort, particles, i_box, ij_components)
        if (abort) return

        call this%visit_walls(overlap, deltas%walls_energies, i_box, ij_components, particles)
        if (overlap) return
        call this%visit_short(overlap, deltas%short_energies, i_box, ij_components, particles)
        if (overlap) return
        call this%visit_field(deltas%field_energies, i_box, particles)
        call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, i_box, &
            ij_components, particles)

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            sum(deltas%short_energies + deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(this%acceptation_probability(i_box, ij_components, &
            delta_energy))
        if (success) call this%update_components(i_box, ij_components, particles)
    end subroutine Abstract_metropolis_algorithm

    !> @note The i_actor <-> j_actor term is ignored.
    pure integer function i_exclude_particle(k_component, ij_components, particles) &
        result(i_exclude)
        integer, intent(in) :: k_component, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        if (k_component == ij_components(1)) then
            i_exclude = particles(1)%i
        else if (k_component == ij_components(2)) then
            i_exclude = particles(2)%i
        else
            i_exclude = 0
        end if
    end function i_exclude_particle

!end implementation Abstract_Box_Particles_Swap

!implementation Box_Particles_Transmutation

    subroutine Transmutation_define_swap(this, abort, particles, i_box, ij_components)
        class(Box_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: particles(:)
        integer, intent(in) :: i_box, ij_components(:)

        if (this%mixture%components(ij_components(1), i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if

        particles(1)%i = random_integer(this%mixture%components(ij_components(1), i_box)%&
            num_particles%get())
        particles(1)%position = this%mixture%components(ij_components(1), i_box)%positions%&
            get(particles(1)%i)
        particles(1)%orientation = this%mixture%components(ij_components(1), i_box)%orientations%&
            get(particles(1)%i)
        particles(1)%dipole_moment = this%mixture%components(ij_components(1), i_box)%&
            dipole_moments%get(particles(1)%i)

        particles(2)%i = this%mixture%components(ij_components(2), i_box)%num_particles%get() + 1
        call this%changes%position_copiers(i_box)%copy(particles(2)%position, &
            particles(1)%position, ij_components)
        call this%changes%orientation_copier%copy(particles(2)%orientation, particles(1)%&
            orientation, ij_components)
        particles(2)%dipole_moment = this%mixture%components(ij_components(2), i_box)%&
            dipole_moments%get_norm() * particles(2)%orientation
    end subroutine Transmutation_define_swap

    !> \[
    !>      P[i \in I \to j \in J] = \min \left( 1, \frac{N_I}{N_J + 1} \frac{\rho_J}{\rho_I}
    !>          e^{-\beta \Delta U_{i \to j}} \frac{a_I^{-N_I}}{a_J^{-N_J}} \right)
    !> \]
    pure real(DP) function Transmutation_acceptation_probability(this, i_box, ij_components, &
        delta_energy) result(probability)
        class(Box_Particles_Transmutation), intent(in) :: this
        integer, intent(in) :: i_box, ij_components(:)
        real(DP), intent(in) :: delta_energy

        associate(component_i => this%mixture%components(ij_components(1), i_box), &
            component_j => this%mixture%components(ij_components(2), i_box))

            probability = &
                real(component_i%num_particles%get(), DP) / &
                real(component_j%num_particles%get() + 1, DP) * &
                component_j%chemical_potential%get_density() / &
                component_i%chemical_potential%get_density() * &
                exp(-delta_energy / this%environment%temperature%get()) * &
                component_i%chemical_potential%get_inv_pow_activity() / &
                component_j%chemical_potential%get_inv_pow_activity()
        end associate

        probability = min(1._DP, probability)
    end function Transmutation_acceptation_probability

    subroutine Transmutation_visit_walls(this, overlap, delta_energies, i_box, ij_components, &
        particles)
        class(Box_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i_partner

        do i_partner = size(particles), 1, -1
            call this%environment%visitable_walls(i_box)%visit(overlap, delta_energies(i_partner), &
                particles(i_partner)%position, this%short_interactions%&
                wall_pairs(ij_components(i_partner))%potential)
            if (overlap) return
        end do
        delta_energies(1) = -delta_energies(1)
    end subroutine Transmutation_visit_walls

    subroutine Transmutation_visit_short(this, overlap, delta_energies, i_box, ij_components, &
        particles)
        class(Box_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:, :)
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: k_component, i_partner, i_exclude

        do i_partner = size(particles), 1, -1
            do k_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
                i_exclude = i_exclude_particle(k_component, ij_components, particles)
                call this%short_interactions%cells(i_box)%visitable_cells(k_component, &
                    ij_components(i_partner))%&
                    visit_energy(overlap, delta_energies(k_component, i_partner), &
                        particles(i_partner), visit_different, i_exclude)
                if (overlap) return
            end do
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
    end subroutine Transmutation_visit_short

    subroutine Transmutation_visit_field(this, delta_energies, i_box, particles)
        class(Box_Particles_Transmutation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        delta_energies = &
            [dipoles_field_visit_remove(this%environment%external_fields(i_box), particles(1)), &
            dipoles_field_visit_add(this%environment%external_fields(i_box), particles(2))]
    end subroutine Transmutation_visit_field

    subroutine Transmutation_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
        ij_components, particles)
        class(Box_Particles_Transmutation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:, :)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: k_component, i_partner, i_exclude

        do i_partner = 1, size(particles)
            do k_component = 1, size(this%dipolar_interactions_dynamic(i_box)%real_components, 1)
                i_exclude = i_exclude_particle(k_component, ij_components, particles)
                call this%dipolar_interactions_dynamic(i_box)%&
                    real_components(k_component, ij_components(i_partner))%&
                    component%visit(delta_energies(k_component, i_partner), particles(i_partner), &
                    visit_different, i_exclude)
            end do
            delta_energies(ij_components(i_partner), i_partner) = &
                delta_energies(ij_components(i_partner), i_partner) - &
                this%dipolar_interactions_dynamic(i_box)%self_components(ij_components(i_partner))%&
                component%meet(particles(i_partner)%dipole_moment)
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
        delta_shared_energy = &
            this%dipolar_interactions_dynamic(i_box)%reci_visitor%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1)) + &
            this%dipolar_interactions_dynamic(i_box)%surf_mixture%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1)%&
                    dipole_moment) - &
            this%dipolar_interactions_dynamic(i_box)%dlc_visitor%&
                visit_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
    end subroutine Transmutation_visit_dipolar

    subroutine Transmutation_update_components(this, i_box, ij_components, particles)
        class(Box_Particles_Transmutation), intent(in) :: this
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) ::  particles(:)

        integer :: k_component

        do k_component = size(this%short_interactions%cells(i_box)%visitable_cells, 2), 1, -1
            call this%short_interactions%cells(i_box)%&
                visitable_cells(ij_components(1), k_component)%remove(particles(1))
        end do

        call this%mixture%total_moments(i_box)%remove(ij_components(1), particles(1)%dipole_moment)
        call this%mixture%components(ij_components(1), i_box)%orientations%remove(particles(1)%i)
        call this%mixture%components(ij_components(1), i_box)%positions%remove(particles(1)%i)
        call this%mixture%components(ij_components(1), i_box)%num_particles%set(this%mixture%&
            components(ij_components(1), i_box)%num_particles%get() - 1)

        call this%mixture%components(ij_components(2), i_box)%num_particles%set(this%mixture%&
            components(ij_components(2), i_box)%num_particles%get() + 1)
        call this%mixture%components(ij_components(2), i_box)%positions%add(particles(2)%position)
        call this%mixture%components(ij_components(2), i_box)%orientations%add(particles(2)%&
            orientation)
        call this%mixture%total_moments(i_box)%add(ij_components(2), particles(2)%dipole_moment)

        do k_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 2)
            call this%short_interactions%cells(i_box)%&
                visitable_cells(ij_components(2), k_component)%add(particles(2))
        end do

        call this%dipolar_interactions_static(i_box)%reci_structure%&
            update_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
        call this%dipolar_interactions_static(i_box)%dlc_structures%&
            update_transmutation(ij_components, particles(2)%dipole_moment, particles(1))
    end subroutine Transmutation_update_components

    subroutine Transmutation_increment_hit(changes, ij_components)
        type(Concrete_Observables_Changes), intent(inout) :: changes
        integer, intent(in) :: ij_components(:)

        changes%transmutations_counters(ij_components(1), ij_components(2))%num_hits = changes%&
            transmutations_counters(ij_components(1), ij_components(2))%num_hits + 1
    end subroutine Transmutation_increment_hit

    subroutine Transmutation_increment_success(changes, ij_components)
        type(Concrete_Observables_Changes), intent(inout) :: changes
        integer, intent(in) :: ij_components(:)

        changes%transmutations_counters(ij_components(1), ij_components(2))%num_successes = &
            changes%transmutations_counters(ij_components(1), ij_components(2))%num_successes + 1
    end subroutine Transmutation_increment_success

!end implementation Box_Particles_Transmutation

!implementation Box_Particles_Switch

    !> @note The change is partially defined because the priority was to share the same interface
    !> with Box_Particles_Transmutation, cf. visit_...() methods for completion.
    subroutine Switch_define_swap(this, abort, particles, i_box, ij_components)
        class(Box_Particles_Switch), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: particles(:)
        integer, intent(in) :: i_box, ij_components(:)

        integer :: i_partner

        if (this%mixture%components(ij_components(1), i_box)%num_particles%get() == 0 .or. &
            this%mixture%components(ij_components(2), i_box)%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if

        do i_partner = 1, size(particles)
            particles(i_partner)%i = random_integer(this%mixture%&
                components(ij_components(i_partner), i_box)%num_particles%get())
            particles(i_partner)%position = this%mixture%&
                components(ij_components(i_partner), i_box)%positions%get(particles(i_partner)%i)
            particles(i_partner)%orientation = this%mixture%&
                components(ij_components(i_partner), i_box)%orientations%get(particles(i_partner)%i)
            particles(i_partner)%dipole_moment = this%mixture%&
                components(ij_components(i_partner), i_box)%dipole_moments%&
                get(particles(i_partner)%i)
        end do
    end subroutine Switch_define_swap

    pure real(DP) function Switch_acceptation_probability(this, i_box, ij_components, delta_energy)&
        result(probability)
        class(Box_Particles_Switch), intent(in) :: this
        integer, intent(in) :: i_box, ij_components(:)
        real(DP), intent(in) :: delta_energy

        probability = min(1._DP, exp(-delta_energy / this%environment%temperature%get()))
    end function Switch_acceptation_probability

    subroutine Switch_visit_walls(this, overlap, delta_energies, i_box, ij_components, particles)
        class(Box_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP) :: energies_new(size(particles)), energies_old(size(particles))
        integer :: i_partner

        do i_partner = 1, size(particles)
            call this%environment%visitable_walls(i_box)%visit(overlap, energies_new(i_partner), &
                particles(size(particles) + 1 - i_partner)%position, this%short_interactions%&
                wall_pairs(ij_components(i_partner))%potential)
            if (overlap) return
        end do
        do i_partner = 1, size(particles)
            call this%environment%visitable_walls(i_box)%visit(overlap, energies_old(i_partner), &
                particles(i_partner)%position, this%short_interactions%&
                wall_pairs(ij_components(i_partner))%potential)
        end do
        delta_energies = energies_new - energies_old
    end subroutine Switch_visit_walls

    subroutine Switch_visit_short(this, overlap, delta_energies, i_box, ij_components, particles)
        class(Box_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:, :)
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP), dimension(size(delta_energies, 1), size(delta_energies, 2)) :: energies_new, &
            energies_old
        type(Concrete_Temporary_Particle) :: swapped(size(particles))
        integer :: k_component, i_partner

        swapped%i = particles%i
        swapped(1)%position = particles(2)%position
        swapped(2)%position = particles(1)%position

        do i_partner = 1, size(swapped)
            do k_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
                call this%short_interactions%cells(i_box)%&
                    visitable_cells(k_component, ij_components(i_partner))%&
                    visit_energy(overlap, energies_new(k_component, i_partner), &
                        swapped(i_partner), visit_different, &
                        i_exclude_particle(k_component, ij_components, swapped))
                if (overlap) return
            end do
        end do
        do i_partner = 1, size(particles)
            do k_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
                call this%short_interactions%cells(i_box)%&
                    visitable_cells(k_component, ij_components(i_partner))%&
                    visit_energy(overlap, energies_old(k_component, i_partner), &
                        particles(i_partner), visit_different, &
                        i_exclude_particle(k_component, ij_components, particles))
            end do
        end do
        delta_energies = energies_new - energies_old
    end subroutine Switch_visit_short

    subroutine Switch_visit_field(this, delta_energies, i_box, particles)
        class(Box_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i_partner

        do i_partner = 1, size(delta_energies)
            delta_energies(i_partner) = dipoles_field_visit_translation(this%environment%&
                external_fields(i_box), particles(size(delta_energies) + 1 - i_partner)%position, &
                particles(i_partner))
        end do
    end subroutine Switch_visit_field

    subroutine Switch_visit_dipolar(this, delta_energies, delta_shared_energy, i_box, &
        ij_components, particles)
        class(Box_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:, :)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP), dimension(size(delta_energies, 1), size(delta_energies, 2)) :: real_energies_new,&
            real_energies_old
        type(Concrete_Temporary_Particle) :: swapped(size(particles))
        integer :: k_component, i_partner

        swapped%i = particles%i
        swapped(1)%position = particles(2)%position
        swapped(2)%position = particles(1)%position
        do i_partner = 1, size(swapped)
            swapped(i_partner)%dipole_moment = particles(i_partner)%dipole_moment
        end do

        do i_partner = 1, size(swapped)
            do k_component = 1, size(this%dipolar_interactions_dynamic(i_box)%real_components, 1)
                call this%dipolar_interactions_dynamic(i_box)%&
                    real_components(k_component, ij_components(i_partner))%component%&
                    visit(real_energies_new(k_component, i_partner), swapped(i_partner), &
                        visit_different, i_exclude_particle(k_component, ij_components, swapped))
                call this%dipolar_interactions_dynamic(i_box)%&
                    real_components(k_component, ij_components(i_partner))%component%&
                    visit(real_energies_old(k_component, i_partner), particles(i_partner), &
                        visit_different, i_exclude_particle(k_component, ij_components, particles))
            end do
        end do
        delta_energies = real_energies_new - real_energies_old
        delta_shared_energy = &
            this%dipolar_interactions_dynamic(i_box)%reci_visitor%&
                visit_switch(ij_components, particles) - &
            this%dipolar_interactions_dynamic(i_box)%dlc_visitor%visit_switch(ij_components, &
                particles)
    end subroutine Switch_visit_dipolar

    subroutine Switch_update_components(this, i_box, ij_components, particles)
        class(Box_Particles_Switch), intent(in) :: this
        integer, intent(in) :: i_box, ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: k_component, i_partner
        real(DP) :: swapped_position(num_dimensions)

        do i_partner = 1, size(particles)
            swapped_position = particles(size(particles) + 1 - i_partner)%position
            call this%mixture%components(ij_components(i_partner), i_box)%positions%&
                set(particles(i_partner)%i, swapped_position)
            do k_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 2)
                call this%short_interactions%cells(i_box)%&
                    visitable_cells(ij_components(i_partner), k_component)%&
                    translate(swapped_position, particles(i_partner))
            end do
        end do

        call this%dipolar_interactions_static(i_box)%reci_structure%&
            update_switch(ij_components, particles)
        call this%dipolar_interactions_static(i_box)%dlc_structures%&
            update_switch(ij_components, particles)
    end subroutine Switch_update_components

    subroutine Switch_increment_hit(changes, ij_components)
        type(Concrete_Observables_Changes), intent(inout) :: changes
        integer, intent(in) :: ij_components(:)

        changes%switches_counters(ij_components(1))%line(ij_components(2))%num_hits = changes%&
            switches_counters(ij_components(1))%line(ij_components(2))%num_hits + 1
    end subroutine Switch_increment_hit

    subroutine Switch_increment_success(changes, ij_components)
        type(Concrete_Observables_Changes), intent(inout) :: changes
        integer, intent(in) :: ij_components(:)

        changes%switches_counters(ij_components(1))%line(ij_components(2))%num_successes = changes%&
            switches_counters(ij_components(1))%line(ij_components(2))%num_successes + 1
    end subroutine Switch_increment_success

!end implementation Box_Particles_Switch

end module classes_box_particles_swap
