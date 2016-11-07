module classes_boxes_particles_switch

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
use procedures_swap_visitors, only: swap_transmutation_visit_short => transmutation_visit_short, &
    swap_transmutation_visit_dipolar => transmutation_visit_dipolar

implicit none

private

    type, extends(Abstract_Generating_Algorithm), public :: Boxes_Particles_Switch
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
        procedure, private :: define_switch => Concrete_define_switch
        procedure, private :: visit_walls => Concrete_visit_walls
        procedure, private :: visit_short => Concrete_visit_short
        procedure, private :: visit_fields => Concrete_visit_fields
        procedure, private :: visit_dipolar => Concrete_visit_dipolar
    end type Boxes_Particles_Switch

contains

    subroutine Concrete_construct(this, environment, mixture, short_interactions, &
        dipolar_interactions_dynamic, dipolar_interactions_static, can_translate, &
        box_couples, component_couples, boxes_selector, components_selectors)
        class(Boxes_Particles_Switch), intent(out) :: this
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
        class(Boxes_Particles_Switch), intent(inout) :: this

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
        class(Boxes_Particles_Switch), intent(inout) :: this

        call selectors_reset(this%components_selectors, this%box_couples, this%component_couples, &
            this%mixture%components, this%can_translate)
    end subroutine Concrete_reset_selectors

    pure integer function Concrete_get_num_choices(this) result(num_choices)
        class(Boxes_Particles_Switch), intent(in) :: this

        integer :: box_i_couple

        num_choices = 0
        do box_i_couple = 1, this%box_couples%get_num()
            num_choices = num_choices + this%components_selectors(box_i_couple)%get_num_choices()
        end do
    end function Concrete_get_num_choices

    subroutine Concrete_try(this, observables)
        class(Boxes_Particles_Switch), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Double_Energies) :: deltas(2)
        integer :: box_i_couple, ij_boxes(size(deltas))
        integer :: i_partner, ij_components(size(deltas))

        box_i_couple = this%boxes_selector%get()
        ij_boxes = this%box_couples%get(box_i_couple)
        ij_components = this%component_couples(box_i_couple)%&
            get(this%components_selectors(box_i_couple)%get())
        observables%switches_counters(maxval(ij_boxes))%line(minval(ij_boxes))%&
            triangle(maxval(ij_components))%line(minval(ij_components))%num_hits = &
            observables%switches_counters(maxval(ij_boxes))%line(minval(ij_boxes))%&
                triangle(maxval(ij_components))%line(minval(ij_components))%num_hits + 1
        do i_partner = 1, size(deltas)
            allocate(deltas(i_partner)%&
                short_energies(size(observables%energies(ij_boxes(i_partner))%short_energies), 2))
            allocate(deltas(i_partner)%&
                dipolar_energies(size(observables%energies(ij_boxes(i_partner))%dipolar_energies), &
                2))
        end do
    end subroutine Concrete_try

    subroutine Concrete_metropolis_algorithm(this, success, deltas, ij_boxes, ij_components)
        class(Boxes_Particles_Switch), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)

        real(DP) :: delta_energies(size(deltas))
        type(Concrete_Double_Particle) :: partners(size(deltas))
        logical :: abort, overlap
        integer :: i_partner

        success = .false.
        call this%define_switch(abort, partners, ij_boxes, ij_components)
        if (abort) return

        call this%visit_walls(overlap, deltas, ij_boxes, ij_components, partners)
        if (overlap) return
        call this%visit_short(overlap, deltas, ij_boxes, ij_components, partners)
        if (overlap) return
        call this%visit_fields(deltas, ij_boxes, partners)
        call this%visit_dipolar(deltas, ij_boxes, ij_components, partners)
    end subroutine Concrete_metropolis_algorithm

    subroutine Concrete_define_switch(this, abort, partners, ij_boxes, ij_components)
        class(Boxes_Particles_Switch), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Double_Particle), intent(inout) :: partners(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)

        integer :: i_partner, i_switched

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
            i_switched = size(partners) + 1 - i_partner
            partners(i_partner)%particles(2)%i = this%mixture%&
                components(ij_components(i_switched), ij_boxes(i_partner))%num_particles%get() + 1
            partners(i_partner)%particles(2)%position = partners(i_partner)%particles(1)%position
            partners(i_partner)%particles(2)%orientation = partners(i_switched)%particles(1)%&
                orientation
            partners(i_partner)%particles(2)%dipole_moment = partners(i_switched)%particles(1)%&
                dipole_moment
        end do
    end subroutine Concrete_define_switch

    pure subroutine Concrete_visit_walls(this, overlap, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call swap_transmutation_visit_short(overlap, deltas(i_partner)%walls_energies, this%&
                environment%visitable_walls(ij_boxes(i_partner)), this%short_interactions%&
                wall_pairs, cshift(ij_components, i_partner-1), partners(i_partner)%particles)
            if (overlap) return
        end do
    end subroutine Concrete_visit_walls

    subroutine Concrete_visit_short(this, overlap, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call swap_transmutation_visit_short(overlap, deltas(i_partner)%short_energies, this%&
                short_interactions%cells(ij_boxes(i_partner)), cshift(ij_components, i_partner-1), &
                partners(i_partner)%particles)
            if (overlap) return
        end do
    end subroutine Concrete_visit_short

    pure subroutine Concrete_visit_fields(this, deltas, ij_boxes, partners)
        class(Boxes_Particles_Switch), intent(in) :: this
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call swap_transmutation_visit_dipolar(deltas(i_partner)%field_energies, this%&
                environment%external_fields(ij_boxes(i_partner)), partners(i_partner)%particles)
        end do
    end subroutine Concrete_visit_fields

    pure subroutine Concrete_visit_dipolar(this, deltas, ij_boxes, ij_components, partners)
        class(Boxes_Particles_Switch), intent(in) :: this
        type(Concrete_Double_Energies), intent(inout) :: deltas(:)
        integer, intent(in) :: ij_boxes(:), ij_components(:)
        type(Concrete_Double_Particle), intent(in) :: partners(:)

        integer :: i_partner

        do i_partner = 1, size(partners)
            call swap_transmutation_visit_dipolar(deltas(i_partner)%dipolar_energies, &
                deltas(i_partner)%dipolar_shared_energy, this%&
                dipolar_interactions_dynamic(ij_boxes(i_partner)), &
                cshift(ij_components, i_partner-1), partners(i_partner)%particles)
        end do
    end subroutine Concrete_visit_dipolar

end module classes_boxes_particles_switch
