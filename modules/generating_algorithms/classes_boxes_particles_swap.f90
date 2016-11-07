module classes_boxes_particles_swap

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use classes_hetero_couples, only: Abstract_Hetero_Couples
use procedures_hetero_couples_factory, only: hetero_couples_destroy => destroy
use classes_tower_sampler, only: Abstract_Tower_Sampler
use procedures_tower_sampler_factory, only: tower_sampler_destroy => destroy
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_particle_wrapper, only: Concrete_Double_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use types_dipolar_interactions_static_wrapper, only: Dipolar_Interactions_Static_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add, &
    dipoles_field_visit_remove => visit_remove
use types_observables_energies, only: Concrete_Double_Energies
use procedures_observables_energies_factory, only: observables_energies_set => set
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_generating_algorithm, only: Abstract_Generating_Algorithm
use procedures_selectors_resetters, only: selectors_reset => reset
use procedures_transmutation_visitors, only: transmutation_visit_short => visit_short, &
    transmutation_visit_dipolar => visit_dipolar
use procedures_metropolis_algorithm, only: metropolis_algorithm
use procedures_transmutation_updaters, only: transmutation_update => update

implicit none

private

    type, extends(Abstract_Generating_Algorithm), public :: Boxes_Particles_Swap
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic(:) => &
            null()
        type(Dipolar_Interactions_Static_Wrapper), pointer :: dipolar_interactions_static(:) => &
            null()
        logical, allocatable :: can_translate(:, :)
        class(Abstract_Hetero_Couples), allocatable :: box_couples, component_couples(:)
        class(Abstract_Tower_Sampler), allocatable :: boxes_selector, components_selectors(:)
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: reset_selectors => Concrete_reset_selectors
        procedure :: get_num_choices => Concrete_get_num_choices
        procedure :: try => Concrete_try
        procedure, private :: metropolis_algorithm => Concrete_metropolis_algorithm
        procedure, private :: define_swap => Concrete_define_swap
        procedure, private :: acceptation_probability => Concrete_acceptation_probability
        procedure, private :: visit_walls => Concrete_visit_walls
        procedure, private :: visit_short => Concrete_visit_short
        procedure, private :: visit_fields => Concrete_visit_fields
        procedure, private :: visit_dipolar => Concrete_visit_dipolar
        procedure, private :: update_components => Concrete_update_components
    end type Boxes_Particles_Swap

contains

    subroutine Concrete_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, can_translate, &
        box_couples, component_couples, boxes_selector, components_selectors)
        class(Boxes_Particles_Swap), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic(:)
        type(Dipolar_Interactions_Static_Wrapper), target, intent(in) :: &
            dipolar_interactions_static(:)
        logical, intent(in) :: can_translate(:, :)
        class(Abstract_Hetero_Couples), intent(in) :: box_couples, component_couples(:)
        class(Abstract_Tower_Sampler), intent(in) :: boxes_selector, components_selectors(:)

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        this%dipolar_interactions_static => dipolar_interactions_static
        allocate(this%can_translate, source=can_translate)
        allocate(this%box_couples, source=box_couples)
        allocate(this%component_couples, source=component_couples)
        allocate(this%boxes_selector, source=boxes_selector)
        allocate(this%components_selectors, source=components_selectors)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Boxes_Particles_Swap), intent(inout) :: this

        call tower_sampler_destroy(this%components_selectors)
        call tower_sampler_destroy(this%boxes_selector)
        call hetero_couples_destroy(this%component_couples)
        call hetero_couples_destroy(this%box_couples)
        if (allocated(this%can_translate)) deallocate(this%can_translate)
        this%dipolar_interactions_static => null()
        this%dipolar_interactions_dynamic => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Concrete_destroy

    subroutine Concrete_reset_selectors(this)
        class(Boxes_Particles_Swap), intent(inout) :: this

        call selectors_reset(this%components_selectors, this%box_couples, this%component_couples, &
            this%mixture%components, this%can_translate)
    end subroutine Concrete_reset_selectors

    pure integer function Concrete_get_num_choices(this) result(num_choices)
        class(Boxes_Particles_Swap), intent(in) :: this

        integer :: box_i_couple

        num_choices = 0
        do box_i_couple = 1, this%box_couples%get_num()
            num_choices = num_choices + this%components_selectors(box_i_couple)%get_num_choices()
        end do
    end function Concrete_get_num_choices

    subroutine Concrete_try(this, observables)
        class(Boxes_Particles_Swap), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Double_Energies) :: deltas(2)
        integer :: box_i_couple, ij_boxes(size(deltas))
        integer :: i_partner, ij_components(size(deltas))

        box_i_couple = this%boxes_selector%get()
        ij_boxes = this%box_couples%get(box_i_couple)
        ij_components = this%component_couples(box_i_couple)%&
            get(this%components_selectors(box_i_couple)%get())
        observables%swaps_counters(ij_components(1), ij_components(2), ij_boxes(1), ij_boxes(2))%&
            num_hits = &
            observables%swaps_counters(ij_components(1), ij_components(2), ij_boxes(1), &
                ij_boxes(2))%num_hits + 1
        do i_partner = 1, size(deltas)
            allocate(deltas(i_partner)%&
                short_energies(size(observables%energies(ij_boxes(i_partner))%short_energies), 2))
            allocate(deltas(i_partner)%&
                dipolar_energies(size(observables%energies(ij_boxes(i_partner))%dipolar_energies), &
                2))
        end do

        call this%metropolis_algorithm(success, deltas, ij_boxes, ij_components)

        if (success) then
            do i_partner = 1, size(ij_boxes)
                observables%nums_particles(ij_components(1), ij_boxes(i_partner)) = this%mixture%&
                    components(ij_components(1), ij_boxes(i_partner))%num_particles%get()
                observables%nums_particles(ij_components(2), ij_boxes(i_partner)) = this%mixture%&
                    components(ij_components(2), ij_boxes(i_partner))%num_particles%get()
                call observables_energies_set(observables%energies(ij_boxes(i_partner)), &
                    deltas(i_partner), cshift(ij_components, i_partner-1))
            end do
            observables%swaps_counters(ij_components(1), ij_components(2), ij_boxes(1), &
                ij_boxes(2))%num_successes = &
                observables%swaps_counters(ij_components(1), ij_components(2), ij_boxes(1), &
                    ij_boxes(2))%num_successes + 1
        end if
    end subroutine Concrete_try

    subroutine Concrete_metropolis_algorithm(this, success, deltas, ij_boxes, ij_components)
        class(Boxes_Particles_Swap), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)

        real(DP) :: delta_energies(size(deltas))
        type(Concrete_Double_Particle) :: partners(size(deltas))
        logical :: abort, overlap
        integer :: i_partner

        success = .false.
        call this%define_swap(abort, partners, ij_boxes, ij_components)
        if (abort) return

        call this%visit_walls(overlap, deltas, ij_boxes, ij_components, partners)
        if (overlap) return
        call this%visit_short(overlap, deltas, ij_boxes, ij_components, partners)
        if (overlap) return
        call this%visit_fields(deltas, ij_boxes, partners)
        call this%visit_dipolar(deltas, ij_boxes, ij_components, partners)

        do i_partner = 1, size(partners)
            delta_energies(i_partner) = &
                sum(deltas(i_partner)%walls_energies + deltas(i_partner)%field_energies) + &
                sum(deltas(i_partner)%short_energies + deltas(i_partner)%dipolar_energies) + &
                deltas(i_partner)%dipolar_shared_energy
        end do
        success = metropolis_algorithm(this%acceptation_probability(ij_boxes, ij_components, &
            sum(delta_energies)))
        if (success) call this%update_components(ij_boxes, ij_components, partners)
    end subroutine Concrete_metropolis_algorithm

    subroutine Concrete_define_swap(this, abort, partners, ij_boxes, ij_components)
        class(Boxes_Particles_Swap), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Double_Particle), intent(inout) :: partners(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)

        integer :: i_partner, i_swapped

        if (this%mixture%components(ij_components(1), ij_boxes(1))%num_particles%get() == 0 .or. &
            this%mixture%components(ij_components(2), ij_boxes(2))%num_particles%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if

        do i_partner = 1, size(partners)
            partners(i_partner)%particles(1)%i = random_integer(this%mixture%&
                components(ij_components(i_partner), ij_boxes(i_partner))%num_particles%get())
            partners(i_partner)%particles(1)%position = this%mixture%&
                components(ij_components(i_partner), ij_boxes(i_partner))%positions%&
                    get(partners(i_partner)%particles(1)%i)
            partners(i_partner)%particles(1)%orientation = this%mixture%&
                components(ij_components(i_partner), ij_boxes(i_partner))%orientations%&
                    get(partners(i_partner)%particles(1)%i)
            partners(i_partner)%particles(1)%dipole_moment = this%mixture%&
                components(ij_components(i_partner), ij_boxes(i_partner))%dipole_moments%&
                    get(partners(i_partner)%particles(1)%i)
        end do

        do i_partner = 1, size(partners)
            i_swapped = size(partners) + 1 - i_partner
            partners(i_partner)%particles(2)%i = this%mixture%&
                components(ij_components(i_swapped), ij_boxes(i_partner))%num_particles%get() + 1
            partners(i_partner)%particles(2)%position = partners(i_partner)%particles(1)%position
            partners(i_partner)%particles(2)%orientation = partners(i_swapped)%particles(1)%&
                orientation
            partners(i_partner)%particles(2)%dipole_moment = partners(i_swapped)%particles(1)%&
                dipole_moment
        end do
    end subroutine Concrete_define_swap

    !> \[
    !>      P[i_I^{\boldsymbol{I}} \to N_I^{\boldsymbol{J}} + 1,
    !>          i_J^{\boldsymbol{J}} \to N_J^{\boldsymbol{I}} + 1] = \min \left(1,
    !>              \frac{N_I^{\boldsymbol{I}}}{N_I^{\boldsymbol{J}} + 1}
    !>              \frac{N_J^{\boldsymbol{J}}}{N_J^{\boldsymbol{I}} + 1}
    !>              e^{-\beta \Delta U_{i_I^{\boldsymbol{I}} \to N_I^{\boldsymbol{J}} + 1,
    !>                  i_J^{\boldsymbol{J}} \to N_J^{\boldsymbol{I}} + 1}} \right)
    !> \]
    pure real(DP) function Concrete_acceptation_probability(this, ij_boxes, ij_components, &
        delta_energy) result(probability)
        class(Boxes_Particles_Swap), intent(in) :: this
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        real(DP), intent(in) :: delta_energy

        probability = &
            real(this%mixture%components(ij_components(1), ij_boxes(1))%num_particles%get(), DP) / &
            real(this%mixture%components(ij_components(1), ij_boxes(2))%num_particles%get()+1, DP)*&
            real(this%mixture%components(ij_components(2), ij_boxes(2))%num_particles%get(), DP) / &
            real(this%mixture%components(ij_components(2), ij_boxes(1))%num_particles%get()+1, DP)*&
            exp(-delta_energy / this%environment%temperature%get())
        probability = min(1._DP, probability)
    end function Concrete_acceptation_probability

    pure subroutine Concrete_visit_walls(this, overlap, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Swap), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call transmutation_visit_short(overlap, deltas(i_partner)%walls_energies, this%&
                environment%visitable_walls(ij_boxes(i_partner)), this%short_interactions%&
                wall_pairs, cshift(ij_components, i_partner-1), partners(i_partner)%particles)
            if (overlap) return
        end do
    end subroutine Concrete_visit_walls

    subroutine Concrete_visit_short(this, overlap, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Swap), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call transmutation_visit_short(overlap, deltas(i_partner)%short_energies, this%&
                short_interactions%cells(ij_boxes(i_partner)), cshift(ij_components, i_partner-1), &
                partners(i_partner)%particles)
            if (overlap) return
        end do
    end subroutine Concrete_visit_short

    pure subroutine Concrete_visit_fields(this, deltas, ij_boxes, partners)
        class(Boxes_Particles_Swap), intent(in) :: this
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call transmutation_visit_dipolar(deltas(i_partner)%field_energies, this%&
                environment%external_fields(ij_boxes(i_partner)), partners(i_partner)%particles)
        end do
    end subroutine Concrete_visit_fields

    pure subroutine Concrete_visit_dipolar(this, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Swap), intent(in) :: this
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call transmutation_visit_dipolar(deltas(i_partner)%dipolar_energies, &
                deltas(i_partner)%dipolar_shared_energy, this%&
                dipolar_interactions_dynamic(ij_boxes(i_partner)), &
                cshift(ij_components, i_partner-1), partners(i_partner)%particles)
        end do
    end subroutine Concrete_visit_dipolar

    subroutine Concrete_update_components(this, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Swap), intent(in) :: this
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call transmutation_update(this%mixture%components(:, ij_boxes(i_partner)), this%&
                mixture%total_moments(ij_boxes(i_partner)), this%short_interactions%&
                cells(ij_boxes(i_partner)), this%dipolar_interactions_static(ij_boxes(i_partner)), &
                cshift(ij_components, i_partner-1), partners(i_partner)%particles)
        end do
    end subroutine Concrete_update_components

end module classes_boxes_particles_swap
