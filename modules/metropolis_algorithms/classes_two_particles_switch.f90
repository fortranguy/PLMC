module classes_two_particles_switch

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_translation => visit_translation
use classes_tower_sampler, only: Abstract_Tower_Sampler
use classes_hetero_couples, only: Abstract_Hetero_Couples
use types_temporary_observables, only: Concrete_Double_Delta_Energies
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use procedures_metropolis_micro, only: update_energies
use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_Two_Particles_Switch
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Wrapper), pointer :: dipolar_interactions => null()
        class(Abstract_Hetero_Couples), allocatable :: couples
        class(Abstract_Tower_Sampler), allocatable :: selector ![i, j] <-> k: convert
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set_selector => Abstract_set_selector
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: test_metropolis => Abstract_test_metropolis
        procedure, private :: define_switch => Abstract_define_switch
        procedure, private :: visit_field => Abstract_visit_field
        procedure, private :: visit_walls => Abstract_visit_walls
        procedure, private :: visit_short => Abstract_visit_short
        procedure, private :: visit_dipolar => Abstract_visit_dipolar
        procedure, private :: update_actors => Abstract_update_actors
    end type Abstract_Two_Particles_Switch

    type, extends(Abstract_Two_Particles_Switch), public :: Concrete_Two_Particles_Switch

    end type Concrete_Two_Particles_Switch

    type, extends(Abstract_Two_Particles_Switch), public :: Null_Two_Particles_Switch
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set_selector => Null_set_selector
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
        procedure, private :: test_metropolis => Null_test_metropolis
        procedure, private :: define_switch => Null_define_switch
        procedure, private :: visit_field => Null_visit_field
        procedure, private :: visit_walls => Null_visit_walls
        procedure, private :: visit_short => Null_visit_short
        procedure, private :: visit_dipolar => Null_visit_dipolar
        procedure, private :: update_actors => Null_update_actors
    end type Null_Two_Particles_Switch

contains

!implementation Abstract_Two_Particles_Switch

    subroutine Abstract_construct(this, environment, components, short_interactions, &
        dipolar_interactions, couples, selector_mold)
        class(Abstract_Two_Particles_Switch), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector_mold

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        this%dipolar_interactions => dipolar_interactions
        allocate(this%couples, source=couples)
        allocate(this%selector, mold=selector_mold)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Two_Particles_Switch), intent(inout) :: this

        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
        if (allocated(this%couples)) then
            call this%couples%destroy()
            deallocate(this%couples)
        end if
        this%dipolar_interactions => null()
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_set_selector(this)
        class(Abstract_Two_Particles_Switch), intent(inout) :: this

        integer :: nums_candidates(this%couples%get_num_indices())
        integer :: i_candidate, ij_couple(2)

        do i_candidate = 1, size(nums_candidates)
            ij_couple = this%couples%get(i_candidate)
            nums_candidates(i_candidate) = minval([this%components(ij_couple(1))%average_number%&
                get(), this%components(ij_couple(2))%average_number%get()])
                !What is the best compromise: minval(), maxval() or average()?
        end do
        call this%selector%construct(nums_candidates)
    end subroutine Abstract_set_selector

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Two_Particles_Switch), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables

        logical :: success
        type(Concrete_Double_Delta_Energies) :: deltas
        integer :: ij_actors(2), i

        ij_actors = this%couples%get(this%selector%get())
        observables%switches_counters(ij_actors(1))%line(ij_actors(2))%num_hits = &
            observables%switches_counters(ij_actors(1))%line(ij_actors(2))%num_hits + 1
        allocate(deltas%short(size(observables%short_energies), 2))
        allocate(deltas%dipolar(size(observables%dipolar_energies), 2))
        call this%test_metropolis(success, deltas, ij_actors)
        if (success) then
            do i = 1, size(ij_actors)
                observables%field_energies(ij_actors(i)) = &
                    observables%field_energies(ij_actors(i)) + deltas%field(i)
                observables%walls_energies(ij_actors(i)) = &
                    observables%walls_energies(ij_actors(i)) + deltas%walls(i)
                call update_energies(observables%short_energies, deltas%short(:, i), ij_actors(i))
                call update_energies(observables%dipolar_energies, deltas%dipolar(:, i), &
                    ij_actors(i))
            end do
            observables%dipolar_mixture_energy = observables%dipolar_mixture_energy + &
                deltas%dipolar_mixture
            observables%switches_counters(ij_actors(1))%line(ij_actors(2))%num_success = &
                observables%switches_counters(ij_actors(1))%line(ij_actors(2))%num_success + 1
        end if
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, deltas, ij_actors)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Delta_Energies), intent(inout) :: deltas
        integer, intent(in) :: ij_actors(:)

        type(Concrete_Temporary_Particle) :: new(2), old(2)
        real(DP) :: delta_energy
        logical :: abort, overlap
        real(DP) :: rand

        success = .false.
        call this%define_switch(abort, new, old, ij_actors)
        if (abort) return

        call this%visit_field(deltas%field, new, old)
        call this%visit_walls(overlap, deltas%walls, ij_actors, new, old)
        if (overlap) return
        call this%visit_short(overlap, deltas%short, ij_actors, new, old)
        if (overlap) return
        call this%visit_dipolar(deltas%dipolar, deltas%dipolar_mixture, ij_actors, new, old)

        delta_energy = sum(deltas%field + deltas%walls) + sum(deltas%short + deltas%dipolar) + &
            deltas%dipolar_mixture
        call random_number(rand)
        if (rand < exp(-delta_energy/this%environment%temperature%get())) then
            call this%update_actors(ij_actors, new, old)
            success = .true.
        end if
    end subroutine Abstract_test_metropolis

    subroutine Abstract_define_switch(this, abort, new, old, ij_actors)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new(:), old(:)
        integer, intent(in) :: ij_actors(:)

        integer :: i

        if (this%components(ij_actors(1))%number%get() == 0 .or. &
            this%components(ij_actors(2))%number%get() == 0) then
            abort = .true.
            return
        else
            abort = .false.
        end if

        do i = 1, size(old)
            old(i)%i = random_integer(this%components(ij_actors(i))%number%get())
            old(i)%position = this%components(ij_actors(i))%positions%get(old(i)%i)
            old(i)%orientation = this%components(ij_actors(i))%orientations%get(old(i)%i)
            old(i)%dipolar_moment = this%components(ij_actors(i))%dipolar_moments%get(old(i)%i)
            new(i)%i = old(i)%i
            new(i)%orientation = old(i)%orientation
            new(i)%dipolar_moment = old(i)%dipolar_moment
        end do

        new(1)%position = old(2)%position
        new(2)%position = old(1)%position
    end subroutine Abstract_define_switch

    subroutine Abstract_visit_field(this, deltas, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        integer :: i

        do i = 1, size(deltas)
            deltas(i) = dipoles_field_visit_translation(this%environment%external_field, new(i)%&
                position, old(i))
        end do
    end subroutine Abstract_visit_field

    subroutine Abstract_visit_walls(this, overlap, deltas, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP) :: energies_new(size(new)), energies_old(size(old))
        integer :: i

        do i = 1, size(new)
            call this%environment%walls%visit(overlap, energies_new(i), new(i)%position, &
                this%short_interactions%wall_pairs(ij_actors(i))%potential)
            if (overlap) return
        end do
        do i = 1, size(old)
            call this%environment%walls%visit(overlap, energies_old(i), old(i)%position, &
                this%short_interactions%wall_pairs(ij_actors(i))%potential)
        end do
        deltas = energies_new - energies_old
    end subroutine Abstract_visit_walls

    subroutine Abstract_visit_short(this, overlap, deltas, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:, :)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP), dimension(size(deltas, 1), size(deltas, 2)) :: energies_new, energies_old
        integer :: i_component, i_exclude, i

        do i = 1, size(new)
            do i_component = 1, size(this%short_interactions%visitable_cells, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, new)
                call this%short_interactions%visitable_cells(i_component, ij_actors(i))%&
                    visit(overlap, energies_new(i_component, i), new(i), i_exclude)
                if (overlap) return
            end do
        end do
        do i = 1, size(old)
            do i_component = 1, size(this%short_interactions%visitable_cells, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, old)
                call this%short_interactions%visitable_cells(i_component, ij_actors(i))%&
                    visit(overlap, energies_old(i_component, i), old(i), i_exclude)
            end do
        end do
        deltas = energies_new - energies_old
    end subroutine Abstract_visit_short

    subroutine Abstract_visit_dipolar(this, deltas, mixture_delta, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: deltas(:, :)
        real(DP), intent(out) :: mixture_delta
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP), dimension(size(deltas, 1), size(deltas, 2)) :: real_energies_new, &
            real_energies_old
        integer :: i_component, i_exclude, i

        do i = 1, size(new)
            do i_component = 1, size(this%dipolar_interactions%real_components, 1)
                !i_actor <-> j_actor: missing
                i_exclude = i_exclude_particle(i_component, ij_actors, new)
                call this%dipolar_interactions%real_components(i_component, ij_actors(i))%&
                    component%visit(real_energies_new(i_component, i), new(i), i_exclude)
                i_exclude = i_exclude_particle(i_component, ij_actors, old)
                call this%dipolar_interactions%real_components(i_component, ij_actors(i))%&
                    component%visit(real_energies_old(i_component, i), old(i), i_exclude)
            end do
        end do
        deltas = real_energies_new - real_energies_old
        mixture_delta = &
            this%dipolar_interactions%reci_visitor%visit_switch(ij_actors, old) - &
            this%dipolar_interactions%dlc_visitor%visit_switch(ij_actors, old)
    end subroutine Abstract_visit_dipolar

    !> Warning: the i_actor <-> j_actor term is ignored.
    pure integer function i_exclude_particle(i_component, ij_actors, particles) result(i_exclude)
        integer, intent(in) :: i_component, ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        if (i_component == ij_actors(1)) then
            i_exclude = particles(1)%i
        else if (i_component == ij_actors(2)) then
            i_exclude = particles(2)%i
        else
            i_exclude = 0
        end if
    end function i_exclude_particle

    subroutine Abstract_update_actors(this, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        integer :: i_component, i

        do i = 1, size(new)
            call this%components(ij_actors(i))%positions%set(new(i)%i, new(i)%position)
            do i_component = 1, size(this%short_interactions%visitable_cells, 2)
                call this%short_interactions%visitable_cells(ij_actors(i), i_component)%&
                    translate(new(i)%position, old(i))
            end do
        end do
        call this%dipolar_interactions%reci_structure%update_switch(ij_actors, old)
        call this%dipolar_interactions%dlc_structures%update_switch(ij_actors, old)
    end subroutine Abstract_update_actors

!end implementation Abstract_Two_Particles_Switch

!implementation Null_Two_Particles_Switch

    subroutine Null_construct(this, environment, components, short_interactions, &
        dipolar_interactions, couples, selector_mold)
        class(Null_Two_Particles_Switch), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        class(Abstract_Hetero_Couples), intent(in) :: couples
        class(Abstract_Tower_Sampler), intent(in) :: selector_mold
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Two_Particles_Switch), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_selector(this)
        class(Null_Two_Particles_Switch), intent(inout) :: this
    end subroutine Null_set_selector

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Two_Particles_Switch), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_Two_Particles_Switch), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

    subroutine Null_test_metropolis(this, success, deltas, ij_actors)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Double_Delta_Energies), intent(inout) :: deltas
        integer, intent(in) :: ij_actors(:)
        success = .false.
        deltas%field = 0._DP; deltas%walls = 0._DP; deltas%short = 0._DP; deltas%dipolar = 0._DP
        deltas%dipolar_mixture = 0._DP
    end subroutine Null_test_metropolis

    subroutine Null_define_switch(this, abort, new, old, ij_actors)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: abort
        type(Concrete_Temporary_Particle), intent(out) :: new(:), old(:)
        integer, intent(in) :: ij_actors(:)
        integer :: i
        abort = .true.
        do i = 1, size(old)
            old(i)%i = 0; new(i)%i = 0
            old(i)%position = 0._DP; new(i)%position = 0._DP
            old(i)%orientation = 0._DP; new(i)%orientation = 0._DP
            old(i)%dipolar_moment = 0._DP; new(i)%dipolar_moment = 0._DP
        end do
    end subroutine Null_define_switch

    subroutine Null_visit_field(this, deltas, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        deltas = 0._DP
    end subroutine Null_visit_field

    subroutine Null_visit_walls(this, overlap, deltas, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        overlap = .false.
        deltas = 0._DP
    end subroutine Null_visit_walls

    subroutine Null_visit_short(this, overlap, deltas, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:, :)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        overlap = .false.
        deltas = 0._DP
    end subroutine Null_visit_short

    subroutine Null_visit_dipolar(this, deltas, mixture_delta, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: deltas(:, :)
        real(DP), intent(out) :: mixture_delta
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        deltas = 0._DP; mixture_delta = 0._DP
    end subroutine Null_visit_dipolar

    subroutine Null_update_actors(this, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
    end subroutine Null_update_actors

!end implementation Null_Two_Particles_Switch

end module classes_two_particles_switch
