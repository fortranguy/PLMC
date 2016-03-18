module class_two_particles_switch

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_random_number, only: random_integer
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_long_interactions_wrapper, only: Long_Interactions_Wrapper
use class_tower_sampler, only: Abstract_Tower_Sampler
use class_hetero_couples, only: Abstract_Hetero_Couples
use module_changes_success, only: Concrete_Switch_Counters
use types_observables_wrapper, only: Observables_Wrapper
use procedures_metropolis_micro, only: update_energies
use class_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_Two_Particles_Switch
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Long_Interactions_Wrapper), pointer :: long_interactions => null()
        class(Abstract_Tower_Sampler), allocatable :: selector ![i, j] <-> k: convert
        class(Abstract_Hetero_Couples), allocatable :: components_couples
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: set_candidates => Abstract_set_candidates
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: try => Abstract_try
        procedure, private :: test_metropolis => Abstract_test_metropolis
        procedure, private :: define_switch => Abstract_define_switch
        procedure, private :: visit_walls => Abstract_visit_walls
        procedure, private :: visit_short => Abstract_visit_short
        procedure, private :: visit_long => Abstract_visit_long
        procedure, private :: update_actors => Abstract_update_actors
    end type Abstract_Two_Particles_Switch

    type, extends(Abstract_Two_Particles_Switch), public :: Concrete_Two_Particles_Switch

    end type Concrete_Two_Particles_Switch

    type, extends(Abstract_Two_Particles_Switch), public :: Null_Two_Particles_Switch
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: set_candidates => Null_set_candidates
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: try => Null_try
        procedure, private :: test_metropolis => Null_test_metropolis
        procedure, private :: define_switch => Null_define_switch
        procedure, private :: visit_walls => Null_visit_walls
        procedure, private :: visit_short => Null_visit_short
        procedure, private :: visit_long => Null_visit_long
        procedure, private :: update_actors => Null_update_actors
    end type Null_Two_Particles_Switch

contains

!implementation Abstract_Two_Particles_Switch

    subroutine Abstract_construct(this, environment, selector, components_couples)
        class(Abstract_Two_Particles_Switch), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Tower_Sampler), intent(in) :: selector
        class(Abstract_Hetero_Couples), intent(in) :: components_couples

        this%environment => environment
        allocate(this%selector, source=selector)
        allocate(this%components_couples, source=components_couples)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Two_Particles_Switch), intent(inout) :: this

        this%long_interactions => null()
        this%short_interactions => null()
        this%components => null()
        if (allocated(this%components_couples)) then
            call this%components_couples%destroy()
            deallocate(this%components_couples)
        end if
        if (allocated(this%selector)) then
            call this%selector%destroy()
            deallocate(this%selector)
        end if
    end subroutine Abstract_destroy

    subroutine Abstract_set_candidates(this, components, short_interactions, long_interactions)
        class(Abstract_Two_Particles_Switch), intent(inout) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), target, intent(in) :: long_interactions

        this%components => components
        this%short_interactions => short_interactions
        this%long_interactions => long_interactions
    end subroutine Abstract_set_candidates

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Two_Particles_Switch), intent(in) :: this

        num_choices = this%selector%get_num_choices()
    end function Abstract_get_num_choices

    subroutine Abstract_try(this, observables)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables

        logical :: success
        real(DP), dimension(2) :: walls_delta, field_energy
        real(DP) :: short_deltas(size(observables%short_energies), 2), &
            long_deltas(size(observables%long_energies), 2)
        real(DP) :: long_mixture_delta
        integer :: ij_actors(2), i

        ij_actors = this%components_couples%get(this%selector%get())
        observables%switches_counters(ij_actors(1))%with_components(ij_actors(2))%num_hits = &
            observables%switches_counters(ij_actors(1))%with_components(ij_actors(2))%num_hits + 1
        call this%test_metropolis(success, walls_delta, short_deltas, long_deltas, &
            long_mixture_delta, ij_actors)
        if (success) then
            do i = 1, size(ij_actors)
                observables%walls_energies(ij_actors(i)) = &
                    observables%walls_energies(ij_actors(i)) + walls_delta(i)
                call update_energies(observables%short_energies, short_deltas(:, i), ij_actors(i))
                call update_energies(observables%long_energies, long_deltas(:, i), ij_actors(i))
            end do
            observables%long_mixture_energy = observables%long_mixture_energy + long_mixture_delta
            observables%switches_counters(ij_actors(1))%with_components(ij_actors(2))%num_success =&
            observables%switches_counters(ij_actors(1))%with_components(ij_actors(2))%num_success +1
        end if
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, walls_delta, short_deltas, long_deltas, &
        long_mixture_delta, ij_actors)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: success
        real(DP), intent(out) :: walls_delta(:)
        real(DP), intent(out) :: short_deltas(:, :), long_deltas(:, :)
        real(DP), intent(out) :: long_mixture_delta
        integer, intent(in) :: ij_actors(:)

        type(Concrete_Temporary_Particle) :: new(2), old(2)
        real(DP) :: energy_delta
        logical :: overlap
        real(DP) :: rand

        call this%define_switch(ij_actors, new, old)

        success = .false.
        call this%visit_walls(overlap, walls_delta, ij_actors, new, old)
        if (overlap) return
        call this%visit_short(overlap, short_deltas, ij_actors, new, old)
        if (overlap) return
        call this%visit_long(long_deltas, long_mixture_delta, ij_actors, new, old)

        energy_delta = sum(walls_delta) + sum(short_deltas + long_deltas) + long_mixture_delta
        call random_number(rand)
        if (rand < exp(-energy_delta/this%environment%temperature%get())) then
            call this%update_actors(ij_actors, new, old)
            success = .true.
        end if
    end subroutine Abstract_test_metropolis

    subroutine Abstract_define_switch(this, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(out) :: new(:), old(:)

        integer :: i

        do i = 1, size(old)
            old(i)%i = random_integer(this%components(ij_actors(i))%positions%get_num())
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

    subroutine Abstract_visit_walls(this, overlap, walls_delta, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: walls_delta(:)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP) :: walls_new(size(new)), walls_old(size(old))
        integer :: i

        do i = 1, size(new)
            call this%environment%walls_potential%visit(overlap, walls_new(i), new(i)%position, &
                this%short_interactions%wall_pairs(ij_actors(i))%pair_potential)
            if (overlap) return
        end do
        do i = 1, size(old)
            call this%environment%walls_potential%visit(overlap, walls_old(i), old(i)%position, &
                this%short_interactions%wall_pairs(ij_actors(i))%pair_potential)
        end do
        walls_delta = walls_new - walls_old
    end subroutine Abstract_visit_walls

    subroutine Abstract_visit_short(this, overlap, short_deltas, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: short_deltas(:, :)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP), dimension(size(short_deltas, 1), size(short_deltas, 2)) :: short_new, short_old
        integer :: i_component, i_exclude, i

        do i = 1, size(new)
            do i_component = 1, size(this%short_interactions%components_cells, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, new)
                call this%short_interactions%components_cells(i_component, ij_actors(i))%&
                    visit(overlap, short_new(i_component, i), new(i), i_exclude)
                if (overlap) return
            end do
        end do
        do i = 1, size(old)
            do i_component = 1, size(this%short_interactions%components_cells, 1)
                i_exclude = i_exclude_particle(i_component, ij_actors, old)
                call this%short_interactions%components_cells(i_component, ij_actors(i))%&
                    visit(overlap, short_old(i_component, i), old(i), i_exclude)
            end do
        end do
        short_deltas = short_new - short_old
    end subroutine Abstract_visit_short

    subroutine Abstract_visit_long(this, long_deltas, long_mixture_delta, ij_actors, new, old)
        class(Abstract_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: long_deltas(:, :)
        real(DP), intent(out) :: long_mixture_delta
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)

        real(DP), dimension(size(long_deltas, 1), size(long_deltas, 2)) :: long_new_real, &
            long_old_real
        integer :: i_component, i_exclude, i

        do i = 1, size(new)
            do i_component = 1, size(this%long_interactions%real_components, 1)
                !i_actor <-> j_actor: missing
                i_exclude = i_exclude_particle(i_component, ij_actors, new)
                call this%long_interactions%real_components(i_component, ij_actors(i))%&
                    real_component%visit(long_new_real(i_component, i), new(i), i_exclude)
                i_exclude = i_exclude_particle(i_component, ij_actors, old)
                call this%long_interactions%real_components(i_component, ij_actors(i))%&
                    real_component%visit(long_old_real(i_component, i), old(i), i_exclude)
            end do
        end do
        long_deltas = long_new_real - long_old_real
        long_mixture_delta = this%long_interactions%reci_visitor%visit_switch(ij_actors, old) - &
            this%long_interactions%dlc_visitor%visit_switch(ij_actors, old)
    end subroutine Abstract_visit_long

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
            do i_component = 1, size(this%short_interactions%components_cells, 1)
                call this%short_interactions%components_cells(ij_actors(i), i_component)%&
                    move(new(i)%position, old(i))
            end do
        end do
        call this%long_interactions%reci_structure%update_switch(ij_actors, old)
        call this%long_interactions%dlc_structures%update_switch(ij_actors, old)
    end subroutine Abstract_update_actors

!end implementation Abstract_Two_Particles_Switch

!implementation Null_Two_Particles_Switch

    subroutine Null_construct(this, environment, selector, components_couples)
        class(Null_Two_Particles_Switch), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Tower_Sampler), intent(in) :: selector
        class(Abstract_Hetero_Couples), intent(in) :: components_couples
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Two_Particles_Switch), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_set_candidates(this, components, short_interactions, long_interactions)
        class(Null_Two_Particles_Switch), intent(inout) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Long_Interactions_Wrapper), target, intent(in) :: long_interactions
    end subroutine Null_set_candidates

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Two_Particles_Switch), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    subroutine Null_try(this, observables)
        class(Null_Two_Particles_Switch), intent(in) :: this
        type(Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

    subroutine Null_test_metropolis(this, success, walls_delta, short_deltas, long_deltas, &
        long_mixture_delta, ij_actors)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: success
        real(DP), intent(out) :: walls_delta(:)
        real(DP), intent(out) :: short_deltas(:, :), long_deltas(:, :)
        real(DP), intent(out) :: long_mixture_delta
        integer, intent(in) :: ij_actors(:)
        success = .false.
        walls_delta = 0._DP; short_deltas = 0._DP; long_deltas = 0._DP; long_mixture_delta = 0._DP
    end subroutine Null_test_metropolis

    subroutine Null_define_switch(this, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(out) :: new(:), old(:)
        integer :: i
        do i = 1, size(old)
            old(i)%i = 0; new(i)%i = 0
            old(i)%position = 0._DP; new(i)%position = 0._DP
            old(i)%orientation = 0._DP; new(i)%orientation = 0._DP
            old(i)%dipolar_moment = 0._DP; new(i)%dipolar_moment = 0._DP
        end do
    end subroutine Null_define_switch

    subroutine Null_visit_walls(this, overlap, walls_delta, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: walls_delta(:)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        overlap = .false.
        walls_delta = 0._DP
    end subroutine Null_visit_walls

    subroutine Null_visit_short(this, overlap, short_deltas, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: short_deltas(:, :)
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        overlap = .false.
        short_deltas = 0._DP
    end subroutine Null_visit_short

    subroutine Null_visit_long(this, long_deltas, long_mixture_delta, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        real(DP), intent(out) :: long_deltas(:, :)
        real(DP), intent(out) :: long_mixture_delta
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
        long_deltas = 0._DP; long_mixture_delta = 0._DP
    end subroutine Null_visit_long

    subroutine Null_update_actors(this, ij_actors, new, old)
        class(Null_Two_Particles_Switch), intent(in) :: this
        integer, intent(in) :: ij_actors(:)
        type(Concrete_Temporary_Particle), intent(in) :: new(:), old(:)
    end subroutine Null_update_actors

!end implementation Null_Two_Particles_Switch

end module class_two_particles_switch
