module classes_two_particles_transmutation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_different
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add, &
    dipoles_field_visit_remove => visit_remove
use classes_tower_sampler, only: Abstract_Tower_Sampler
use classes_hetero_couples, only: Abstract_Hetero_Couples
use types_changes_wrapper, only: Changes_Wrapper
use types_observables_energies, only: Concrete_Double_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_metropolis_algorithm, only: metropolis_algorithm

implicit none

private

    type, extends(Abstract_Generating_Algorithm), abstract, public :: &
        Abstract_Two_Particles_Transmutation
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static => null()
        type(Changes_Wrapper), pointer :: changes => null()
        logical, allocatable :: can_exchange(:)
        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector ![i, j] <-> k: convert
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set_selector => Abstract_set_selector
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: metropolis_algorithm => Abstract_metropolis_algorithm
        procedure, private :: define_transmutation => Abstract_define_transmutation
        procedure, private :: acceptation_probability => Abstract_acceptation_probability
        procedure, private :: visit_walls => Abstract_visit_walls
        procedure, private :: visit_short => Abstract_visit_short
        procedure, private :: visit_field => Abstract_visit_field
        procedure, private :: visit_dipolar => Abstract_visit_dipolar
        procedure, private :: update_actors => Abstract_update_actors
    end type Abstract_Two_Particles_Transmutation

    type, extends(Abstract_Two_Particles_Transmutation), public :: &
        Concrete_Two_Particles_Transmutation

    end type Concrete_Two_Particles_Transmutation

    type, extends(Abstract_Two_Particles_Transmutation), public :: Null_Two_Particles_Transmutation
    contains
    procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set_selector => Null_set_selector
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
    end type Null_Two_Particles_Transmutation

contains

!implementation Abstract_Two_Particles_Transmutation

    subroutine Abstract_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes, can_exchange, couples, &
        selector_mold)
        class(Abstract_Two_Particles_Transmutation), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static
        type(Changes_Wrapper), target, intent(in) :: changes
        logical, intent(in) :: can_exchange(:)
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector_mold

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%changes => changes
        allocate(this%can_exchange, source=can_exchange)
        allocate(this%couples, source=couples)
        allocate(this%selector, mold=selector_mold)
        !this%selector: delayed construction in [[Abstract_set_selector]]
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Two_Particles_Transmutation), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        if (allocated(this%couples)) then
            call this%couples%destroy()
            deallocate(this%couples)
        end if
        !this%can_exchange: early deallocation in [[Abstract_set_selector]]
        this%changes => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set_selector(this)
        class(Abstract_Two_Particles_Transmutation), intent(inout) :: this

        integer :: nums_candidates(this%couples%get_num_indices())
        integer :: i_candidate, ij_couple(2)

        do i_candidate = 1, size(nums_candidates)
            ij_couple = this%couples%get(i_candidate)
            nums_candidates(i_candidate) = &
                merge(minval([this%mixture%components(ij_couple(1))%average_number%get(), &
                    this%mixture%components(ij_couple(2))%average_number%get()]), &
                0, this%can_exchange(ij_couple(1)) .and. this%can_exchange(ij_couple(2)))
                !What is the best compromise: minval(), maxval() or average()?
        end do
        if (allocated(this%can_exchange)) deallocate(this%can_exchange)
        call this%selector%construct(nums_candidates)
    end subroutine Abstract_set_selector

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Double_Energies) :: deltas
        integer :: ij_actors(2)

        ij_actors = this%couples%get(this%selector%get())
        observables%transmutations_counters(ij_actors(1), ij_actors(2))%num_hits = &
            observables%transmutations_counters(ij_actors(1), ij_actors(2))%num_hits + 1
        allocate(deltas%short_energies(size(observables%energies%short_energies), 2))
        allocate(deltas%dipolar_energies(size(observables%energies%dipolar_energies), 2))
        call this%metropolis_algorithm(success, deltas, ij_actors)
        if (success) then
            call observables_energies_set(observables%energies, deltas, ij_actors)
            observables%transmutations_counters(ij_actors(1), ij_actors(2))%num_successes = &
                observables%transmutations_counters(ij_actors(1), ij_actors(2))%num_successes + 1
        end if
    end subroutine Abstract_try

    subroutine Abstract_metropolis_algorithm(this, success, deltas, ij_actors)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Energies), intent(inout) :: deltas
        integer, intent(in) :: ij_actors(:)

        type(Concrete_Temporary_Particle) :: particles(2)
        real(DP) :: delta_energy
        logical :: abort, overlap

        success = .false.
        call this%define_transmutation(abort, particles, ij_actors)
        if (abort) return

        call this%visit_walls(overlap, deltas%walls_energies, ij_actors, particles)
        if (overlap) return
        call this%visit_short(overlap, deltas%short_energies, ij_actors, particles)
        if (overlap) return
        call this%visit_field(deltas%field_energies, particles)
        call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, ij_actors, &
            particles)

        delta_energy = sum(deltas%walls_energies + deltas%field_energies) + &
            sum(deltas%short_energies + deltas%dipolar_energies) + deltas%dipolar_shared_energy
        success = metropolis_algorithm(this%acceptation_probability(ij_actors, delta_energy))
        if (success) call this%update_actors(ij_actors, particles)
    end subroutine Abstract_metropolis_algorithm

    subroutine Abstract_define_transmutation(this, abort, particles, ij_actors)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: particles(:)
        integer, intent(in) :: ij_actors(:)

        if (this%mixture%components(ij_actors(1))%number%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if
        particles(1)%i = random_integer(this%mixture%components(ij_actors(1))%number%get())
        particles(1)%position = this%mixture%components(ij_actors(1))%positions%get(particles(1)%i)
        particles(1)%orientation = this%mixture%components(ij_actors(1))%orientations%&
            get(particles(1)%i)
        particles(1)%dipole_moment = this%mixture%components(ij_actors(1))%dipole_moments%&
            get(particles(1)%i)
        particles(2)%i = this%mixture%components(ij_actors(2))%number%get() + 1
        call this%changes%position_copier%copy(particles(2)%position, particles(1)%position, &
            ij_actors)
        call this%changes%orientation_copier%copy(particles(2)%orientation, particles(1)%&
            orientation, ij_actors)
        particles(2)%dipole_moment = this%mixture%components(ij_actors(2))%dipole_moments%&
            get_norm() * particles(2)%orientation
    end subroutine Abstract_define_transmutation

    !> \[
    !>      P[i \in I \to j \in J] = min \left( 1, \frac{N_I}{N_J + 1} \frac{\rho_J}{\rho_I}
    !>          e^{-\beta \Delta U_{i \to j}} \frac{a_I^{-N_I}}{a_J^{-N_J}} \right)
    !> \]
    pure real(DP) function Abstract_acceptation_probability(this, ij_actors, delta_energy) &
        result(probability)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        real(DP), intent(in) :: delta_energy

        associate(temperature => this%environment%temperature%get(), &
            component_i => this%mixture%components(ij_actors(1)), &
            component_j => this%mixture%components(ij_actors(2)))
            probability = &
                real(component_i%number%get(), DP) / real(component_j%number%get() + 1, DP) * &
                component_j%chemical_potential%get_density() / component_i%chemical_potential%&
                    get_density() * &
                exp(-delta_energy/temperature) * component_i%chemical_potential%&
                get_inv_pow_activity() / component_j%chemical_potential%get_inv_pow_activity()
        end associate
        probability = min(1._DP, probability)
    end function Abstract_acceptation_probability

    subroutine Abstract_visit_walls(this, overlap, delta_energies, ij_actors, particles)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i

        do i = size(particles), 1, -1
            call this%environment%visitable_walls%visit(overlap, delta_energies(i), particles(i)%&
                position, this%short_interactions%wall_pairs(ij_actors(i))%potential)
            if (overlap) return
        end do
        delta_energies(1) = -delta_energies(1)
    end subroutine Abstract_visit_walls

    subroutine Abstract_visit_short(this, overlap, delta_energies, ij_actors, particles)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta_energies(:, :)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i_component, i_exclude, i

        do i = size(particles), 1, -1
            do i_component = 1, size(this%short_interactions%visitable_cells, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, particles, i)
                call this%short_interactions%visitable_cells(i_component, ij_actors(i))%&
                    visit_energy(overlap, delta_energies(i_component, i), particles(i), &
                    visit_different, i_exclude)
                if (overlap) return
            end do
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
    end subroutine Abstract_visit_short

    subroutine Abstract_visit_field(this, delta_energies, particles)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        delta_energies = &
            [dipoles_field_visit_remove(this%environment%external_field, particles(1)), &
            dipoles_field_visit_add(this%environment%external_field, particles(2))]
    end subroutine Abstract_visit_field

    subroutine Abstract_visit_dipolar(this, delta_energies, delta_shared_energy, ij_actors, &
        particles)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        real(DP), intent(out) :: delta_energies(:, :)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i_component, i_exclude, i

        do i = 1, size(particles)
            do i_component = 1, size(this%dipolar_interactions_dynamic%real_components, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, particles, i)
                call this%dipolar_interactions_dynamic%real_components(i_component, ij_actors(i))%&
                    component%visit(delta_energies(i_component, i), particles(i), visit_different, &
                    i_exclude)
            end do
            delta_energies(ij_actors(i), i) = delta_energies(ij_actors(i), i) - &
                this%dipolar_interactions_dynamic%self_components(ij_actors(i))%component%&
                meet(particles(i)%dipole_moment)
        end do
        delta_energies(:, 1) = -delta_energies(:, 1)
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitor%&
                visit_transmutation(ij_actors, particles(2)%dipole_moment, particles(1)) + &
            this%dipolar_interactions_dynamic%surf_mixture%&
                visit_transmutation(ij_actors, particles(2)%dipole_moment, &
                    particles(1)%dipole_moment) - &
            this%dipolar_interactions_dynamic%dlc_visitor%&
                visit_transmutation(ij_actors, particles(2)%dipole_moment, particles(1))
    end subroutine Abstract_visit_dipolar

    pure integer function i_exclude_particle(i_component, ij_actors, particles, i) result(i_exclude)
        integer, intent(in) :: i_component, ij_actors(:), i
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        if (i==1) then
            i_exclude = merge(particles(1)%i, 0, i_component == ij_actors(1))
        else if (i==2) then
            if (i_component == ij_actors(1)) then
                i_exclude = particles(1)%i
            else if (i_component == ij_actors(2)) then
                i_exclude = particles(2)%i
            else
                i_exclude = 0
            end if
        !else: error
        end if
    end function i_exclude_particle

    subroutine Abstract_update_actors(this, ij_actors, particles)
        class(Abstract_Two_Particles_Transmutation), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) ::  particles(:)

        integer :: i_component

        do i_component = size(this%short_interactions%visitable_cells, 2), 1, -1
            call this%short_interactions%visitable_cells(ij_actors(1), i_component)%&
                remove(particles(1))
        end do
        call this%mixture%total_moment%remove(ij_actors(1), particles(1)%dipole_moment)
        call this%mixture%components(ij_actors(1))%orientations%remove(particles(1)%i)
        call this%mixture%components(ij_actors(1))%positions%remove(particles(1)%i)
        call this%mixture%components(ij_actors(1))%number%set(this%mixture%&
            components(ij_actors(1))%number%get() - 1)
        call this%mixture%components(ij_actors(2))%number%set(this%mixture%&
            components(ij_actors(2))%number%get() + 1)
        call this%mixture%components(ij_actors(2))%positions%add(particles(2)%position)
        call this%mixture%components(ij_actors(2))%orientations%add(particles(2)%orientation)
        call this%mixture%total_moment%add(ij_actors(2), particles(2)%dipole_moment)
        do i_component = 1, size(this%short_interactions%visitable_cells, 2)
            call this%short_interactions%visitable_cells(ij_actors(2), i_component)%&
                add(particles(2))
        end do

        call this%dipolar_interactions_static%reci_structure%&
            update_transmutation(ij_actors, particles(2)%dipole_moment, particles(1))
        call this%dipolar_interactions_static%dlc_structures%&
            update_transmutation(ij_actors, particles(2)%dipole_moment, particles(1))
    end subroutine Abstract_update_actors

!end implementation Abstract_Two_Particles_Transmutation

!implementation Null_Two_Particles_Transmutation

    subroutine Null_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, changes, can_exchange, couples, &
        selector_mold)
        class(Null_Two_Particles_Transmutation), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: dipolar_interactions_static
        type(Changes_Wrapper), target, intent(in) :: changes
        logical, intent(in) :: can_exchange(:)
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector_mold
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Two_Particles_Transmutation), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_selector(this)
        class(Null_Two_Particles_Transmutation), intent(inout) :: this
    end subroutine Null_set_selector

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Two_Particles_Transmutation), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_Two_Particles_Transmutation), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Two_Particles_Transmutation

end module classes_two_particles_transmutation
