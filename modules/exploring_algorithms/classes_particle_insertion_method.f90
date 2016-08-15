module classes_particle_insertion_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_environment_wrapper, only: Environment_Wrapper
use classes_component_number, only: Abstract_Component_Number
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_condition_different => different
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add
use classes_random_coordinates, only: Abstract_Random_Coordinates
use types_observables_energies, only: Concrete_Single_Energies
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Particle_Insertion_Method
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Wrapper), pointer :: dipolar_interactions => null()
        class(Abstract_Random_Coordinates), pointer :: random_position => null(), &
            random_orientation => null()
        class(Abstract_Component_Number), allocatable :: numbers(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure, private :: visit_field => Abstract_visit_field
        procedure, private :: visit_walls => Abstract_visit_walls
        procedure, private :: visit_short => Abstract_visit_short
        procedure, private :: visit_dipolar => Abstract_visit_dipolar
    end type Abstract_Particle_Insertion_Method

    type, extends(Abstract_Particle_Insertion_Method), public :: Concrete_Particle_Insertion_Method

    end type Concrete_Particle_Insertion_Method

    type, extends(Abstract_Particle_Insertion_Method), public :: Null_Particle_Insertion_Method
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: try => Null_try
    end type Null_Particle_Insertion_Method

contains

!implementation Abstract_Particle_Insertion_Method

    subroutine Abstract_construct(this, environment, components, short_interactions, &
        dipolar_interactions, numbers, random_position, random_orientation)
        class(Abstract_Particle_Insertion_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        class(Abstract_Component_Number), intent(in) :: numbers(:)
        class(Abstract_Random_Coordinates), target, intent(in) :: random_position, &
            random_orientation

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
        this%dipolar_interactions => dipolar_interactions
        allocate(this%numbers, source=numbers)
        this%random_position => random_position
        this%random_orientation => random_orientation
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Particle_Insertion_Method), intent(inout) :: this

        this%random_orientation => null()
        this%random_position => null()
        if (allocated(this%numbers)) deallocate(this%numbers)
        this%dipolar_interactions => null()
        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        real(DP) :: inv_pow_activity_sum, delta_energy
        type(Concrete_Single_Energies) :: deltas
        type(Concrete_Temporary_Particle) :: test
        integer :: i_component, i_particle
        logical :: overlap

        allocate(deltas%short(size(observables%inv_pow_activities)))
        allocate(deltas%dipolar(size(deltas%short)))
        observables%inv_pow_activities = 0._DP
        test%i = 0
        do i_component = 1, size(this%numbers)
            inv_pow_activity_sum = 0._DP
            do i_particle = 1, this%numbers(i_component)%get()

                observables%insertion_counters(i_component)%num_hits = observables%&
                    insertion_counters(i_component)%num_hits + 1
                test%position = this%random_position%get(i_component)
                test%orientation = this%random_orientation%get(i_component)
                test%dipole_moment = this%components(i_component)%dipole_moments%get_norm() * &
                    test%orientation

                call this%visit_field(deltas%field, test)
                call this%visit_walls(overlap, deltas%walls, i_component, test)
                if (overlap) cycle
                call this%visit_short(overlap, deltas%short, i_component, test)
                if (overlap) cycle
                call this%visit_dipolar(deltas%dipolar, deltas%dipolar_mixture, i_component, test)

                delta_energy = deltas%field + deltas%walls + sum(deltas%short + deltas%dipolar) + &
                    deltas%dipolar_mixture
                inv_pow_activity_sum = inv_pow_activity_sum + exp(-delta_energy/this%environment%&
                    temperature%get())
                observables%insertion_counters(i_component)%num_successes = observables%&
                    insertion_counters(i_component)%num_successes + 1
            end do
            observables%inv_pow_activities(i_component) = inv_pow_activity_sum / real(this%&
                numbers(i_component)%get())
        end do
    end subroutine Abstract_try

    subroutine Abstract_visit_field(this, delta, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        real(DP), intent(out) :: delta
        type(Concrete_Temporary_Particle), intent(in) :: test

        delta = dipoles_field_visit_add(this%environment%external_field, test)
    end subroutine Abstract_visit_field

    subroutine Abstract_visit_walls(this, overlap, delta, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        call this%environment%visitable_walls%visit(overlap, delta, test%position, this%&
            short_interactions%wall_pairs(i_component)%potential)
    end subroutine Abstract_visit_walls

    subroutine Abstract_visit_short(this, overlap, deltas, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        integer :: j_component, i_exclude

        do j_component = 1, size(this%short_interactions%visitable_cells, 1)
            i_exclude = merge(test%i, 0, j_component == i_component)
            call this%short_interactions%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, deltas(j_component), test, visit_condition_different, &
                i_exclude)
                if (overlap) return
        end do
    end subroutine Abstract_visit_short

    subroutine Abstract_visit_dipolar(this, deltas, mixture_delta, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        real(DP), intent(out) :: mixture_delta
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions%real_components, 1)
            i_exclude = merge(test%i, 0, j_component == i_component)
            call this%dipolar_interactions%real_components(j_component, i_component)%component%&
                visit(deltas(j_component), test, i_exclude)
        end do
        deltas(i_component) = deltas(i_component) - this%dipolar_interactions%&
            self_components(i_component)%component%meet(test%dipole_moment)
        mixture_delta = this%dipolar_interactions%reci_visitor%visit_add(i_component, test) + &
            this%dipolar_interactions%surf_mixture%visit_add(i_component, test%dipole_moment) - &
            this%dipolar_interactions%dlc_visitor%visit_add(i_component, test)
    end subroutine Abstract_visit_dipolar

!end implementation Abstract_Particle_Insertion_Method

!implementation Null_Particle_Insertion_Method

    subroutine Null_construct(this, environment, components, short_interactions, &
        dipolar_interactions, numbers, random_position, random_orientation)
        class(Null_Particle_Insertion_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions
        class(Abstract_Component_Number), intent(in) :: numbers(:)
        class(Abstract_Random_Coordinates), target, intent(in) :: random_position, &
            random_orientation
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Particle_Insertion_Method), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, observables)
        class(Null_Particle_Insertion_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Particle_Insertion_Method

end module classes_particle_insertion_method
