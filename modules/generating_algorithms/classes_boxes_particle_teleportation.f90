module classes_boxes_particle_teleportation

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add, &
    dipoles_field_visit_remove => visit_remove
use classes_random_coordinates, only: Abstract_Random_Coordinates
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

    type, extends(Abstract_Generating_Algorithm), public :: Boxes_Particle_Teleportation
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic(:) => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        class(Abstract_Random_Coordinates), pointer :: random_positions(:) => null()
        logical, allocatable :: can_translate(:, :)
        class(Abstract_Hetero_Couples), allocatable :: box_couples
        class(Abstract_Tower_Sampler), allocatable :: boxes_selector, component_selectors(:)
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: reset_selectors => Concrete_reset_selectors
        procedure :: get_num_choices => Concrete_get_num_choices
        procedure :: try => Concrete_try
        procedure, private :: metropolis_algorithm => Concrete_metropolis_algorithm
        procedure, private :: define_teleportation => Concrete_define_teleportation
        procedure, private :: acceptation_probability => Concrete_acceptation_probability
        procedure, private :: visit_walls => Concrete_visit_walls
        procedure, private :: visit_short => Concrete_visit_short
        procedure, private :: visit_field => Concrete_visit_field
        procedure, private :: visit_dipolar => Concrete_visit_dipolar
        procedure, private :: update_components => Concrete_update_components
    end type Boxes_Particle_Teleportation

contains

    subroutine Concrete_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, random_positions, can_translate,&
        box_couples, boxes_selector, component_selectors)
        class(Boxes_Particle_Teleportation), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic(:)
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: &
            dipolar_interactions_static(:)
        class(Abstract_Random_Coordinates), target, intent(in) :: random_positions(:)
        logical, intent(in) :: can_translate(:, :)
        class(Abstract_Hetero_Couples), intent(in) :: box_couples
        class(Abstract_Tower_Sampler), intent(in) :: boxes_selector, component_selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        this%random_positions => random_positions
        allocate(this%can_translate, source=can_translate)
        allocate(this%box_couples, source=box_couples)
        allocate(this%boxes_selector, source=boxes_selector)
        allocate(this%component_selectors, source=component_selectors)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Boxes_Particle_Teleportation), intent(inout) :: this

        call tower_sampler_destroy(this%component_selectors)
        call tower_sampler_destroy(this%boxes_selector)
        call hetero_couples_destroy(this%box_couples)
        if (allocated(this%can_translate)) deallocate(this%can_translate)
        this%random_positions => null()
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Concrete_destroy

    subroutine Concrete_reset_selectors(this)
        class(Boxes_Particle_Teleportation), intent(inout) :: this

        call selectors_reset(this%component_selectors, this%mixture%components, this%can_translate)
    end subroutine Concrete_reset_selectors

    !> @todo Is num_choices enough?
    pure integer function Concrete_get_num_choices(this) result(num_choices)
        class(Boxes_Particle_Teleportation), intent(in) :: this

        integer :: i_box

        num_choices = 0
        do i_box = 1, size(this%component_selectors)
            num_choices = num_choices + this%component_selectors(i_box)%get_num_choices()
        end do
    end function Concrete_get_num_choices

    subroutine Concrete_try(this, observables)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Single_Energies) :: deltas(2)
        integer :: i_partner, ij_boxes(size(deltas)), i_component

        ij_boxes = this%box_couples%get(this%boxes_selector%get())
        i_component = this%component_selectors(ij_boxes(1))%get()
        observables%teleportations_counters(i_component, ij_boxes(1), ij_boxes(2))%num_hits = &
            observables%teleportations_counters(i_component, ij_boxes(1), ij_boxes(2))%num_hits + 1
        do i_partner = 1, size(deltas)
            allocate(deltas(i_partner)%&
                short_energies(size(observables%energies(ij_boxes(i_partner))%short_energies)))
            allocate(deltas(i_partner)%&
                dipolar_energies(size(observables%energies(ij_boxes(i_partner))%dipolar_energies)))
        end do
        call this%metropolis_algorithm(success, deltas, ij_boxes, i_component)
        if (success) then
            do i_partner = 1, size(deltas)
                observables%nums_particles(i_component, ij_boxes(i_partner)) = this%mixture%&
                    components(i_component, ij_boxes(i_partner))%num_particles%get()
                call observables_energies_set(observables%energies(ij_boxes(i_partner)), &
                    deltas(i_partner), i_component)
            end do
            observables%teleportations_counters(i_component, ij_boxes(1), ij_boxes(2))%&
                num_successes = &
                observables%teleportations_counters(i_component, ij_boxes(1), ij_boxes(2))%&
                    num_successes + 1
        end if
    end subroutine Concrete_try

    subroutine Concrete_metropolis_algorithm(this, success, deltas, ij_boxes, i_component)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), i_component

        real(DP) :: delta_energy
        type(Concrete_Temporary_Particle) :: particles(2)
        logical :: abort, overlap
        integer :: i_partner

        success = .false.
        call this%define_teleportation(abort, particles, ij_boxes, i_component)
        if (abort) return

        call this%visit_walls(overlap, deltas, ij_boxes, i_component, particles)
        if (overlap) return
        call this%visit_short(overlap, deltas, ij_boxes, i_component, particles)
        if (overlap) return
        call this%visit_field(deltas, ij_boxes, particles)
        call this%visit_dipolar(deltas, ij_boxes, i_component, particles)

        delta_energy = sum(deltas%walls_energy + deltas%field_energy) + &
            sum(deltas%dipolar_shared_energy)
        do i_partner = 1, size(deltas)
            delta_energy = delta_energy + &
                sum(deltas(i_partner)%short_energies + deltas(i_partner)%dipolar_energies)
        end do
        success = metropolis_algorithm(this%acceptation_probability(ij_boxes, i_component, &
            delta_energy))
        if (success) call this%update_components(ij_boxes, i_component, particles)
    end subroutine Concrete_metropolis_algorithm

    subroutine Concrete_define_teleportation(this, abort, particles, ij_boxes, i_component)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: particles(:)
        integer, intent(in) :: ij_boxes(:), i_component

        if (this%mixture%components(i_component, ij_boxes(1))%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if

        particles(1)%i = random_integer(this%mixture%components(i_component, ij_boxes(1))%&
            num_particles%get())
        particles(1)%position = this%mixture%components(i_component, ij_boxes(1))%positions%&
            get(particles(1)%i)
        particles(1)%orientation = this%mixture%components(i_component, ij_boxes(1))%orientations%&
            get(particles(1)%i)
        particles(1)%dipole_moment = this%mixture%components(i_component, ij_boxes(1))%&
            dipole_moments%get(particles(1)%i)

        particles(2)%i = this%mixture%components(i_component, ij_boxes(2))%num_particles%get() + 1
        particles(2)%position = this%random_positions(ij_boxes(2))%get(i_component)
        particles(2)%orientation = particles(1)%orientation
        particles(2)%dipole_moment = particles(1)%dipole_moment
    end subroutine Concrete_define_teleportation

    !> \[
    !>      P \left[ i_I^{\boldsymbol{I}} \to N_I^{\boldsymbol{J}} + 1 \right] = \min \left(1,
    !>          \frac{N_I^{\boldsymbol{I}}}{N_I^{\boldsymbol{J}} + 1}
    !>          \frac{V^{\boldsymbol{J}}}{V^{\boldsymbol{I}}}
    !>          e^{-\beta \Delta U_{i_I^{\boldsymbol{I}} \to N_I^{\boldsymbol{J}} + 1}} \right)
    !> \]
    pure real(DP) function Concrete_acceptation_probability(this, ij_boxes, i_component, &
        delta_energy) result(probability)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        integer, intent(in) :: ij_boxes(:), i_component
        real(DP), intent(in) :: delta_energy

        probability = &
            real(this%mixture%components(i_component, ij_boxes(1))%num_particles%get(), DP) / &
            real(this%mixture%components(i_component, ij_boxes(2))%num_particles%get() + 1, DP) * &
            product(this%environment%accessible_domains(ij_boxes(2))%get_size()) / &
            product(this%environment%accessible_domains(ij_boxes(1))%get_size()) * &
            exp(-delta_energy / this%environment%temperature%get())

        probability = min(1._DP, probability)
    end function Concrete_acceptation_probability

    subroutine Concrete_visit_walls(this, overlap, deltas, ij_boxes, i_component, particles)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Single_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), i_component
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        integer :: i_partner

        do i_partner = size(ij_boxes), 1, -1
            call this%environment%visitable_walls(ij_boxes(i_partner))%&
                visit(overlap, deltas(i_partner)%walls_energy, particles(i_partner)%position, this%&
                    short_interactions%wall_pairs(i_component)%potential)
            if (overlap) return
        end do
        deltas(1)%walls_energy = -deltas(1)%walls_energy
    end subroutine Concrete_visit_walls

    subroutine Concrete_visit_short(this, overlap, deltas, ij_boxes, i_component, particles)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Single_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), i_component
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        call exchange_visit_add(overlap, deltas(2)%short_energies, i_component, particles(2), this%&
            short_interactions%cells(ij_boxes(2)))
        if (overlap) return
        call exchange_visit_remove(overlap, deltas(1)%short_energies, i_component, particles(1), &
            this%short_interactions%cells(ij_boxes(1)))
    end subroutine Concrete_visit_short

    subroutine Concrete_visit_field(this, deltas, ij_boxes, particles)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        type(Concrete_Single_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        deltas(1)%field_energy = dipoles_field_visit_remove(this%environment%&
            external_fields(ij_boxes(1)), particles(1))
        deltas(2)%field_energy = dipoles_field_visit_add(this%environment%&
            external_fields(ij_boxes(2)), particles(2))
    end subroutine Concrete_visit_field

    subroutine Concrete_visit_dipolar(this, deltas, ij_boxes, i_component, particles)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        type(Concrete_Single_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), i_component
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        call exchange_visit_add(deltas(2)%dipolar_energies, deltas(2)%dipolar_shared_energy, &
            i_component, particles(2), this%dipolar_interactions_dynamic(ij_boxes(2)))
        call exchange_visit_remove(deltas(1)%dipolar_energies, deltas(1)%dipolar_shared_energy, &
            i_component, particles(1), this%dipolar_interactions_dynamic(ij_boxes(1)))
    end subroutine Concrete_visit_dipolar

    subroutine Concrete_update_components(this, ij_boxes, i_component, particles)
        class(Boxes_Particle_Teleportation), intent(in) :: this
        integer, intent(in) :: ij_boxes(:), i_component
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        call exchange_update_remove(this%mixture%components(i_component, ij_boxes(1)), &
            this%mixture%total_moments(ij_boxes(1)), this%short_interactions%cells(ij_boxes(1)), &
            this%dipolar_interactions_static(ij_boxes(1)), i_component, particles(1))

        call exchange_update_add(this%mixture%components(i_component, ij_boxes(2)), &
            this%mixture%total_moments(ij_boxes(2)), this%short_interactions%cells(ij_boxes(2)), &
            this%dipolar_interactions_static(ij_boxes(2)), i_component, particles(2))
    end subroutine Concrete_update_components

end module classes_boxes_particle_teleportation
