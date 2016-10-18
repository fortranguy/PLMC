module classes_particle_insertion_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_environment_wrapper, only: Environment_Wrapper
use classes_num_particles, only: Abstract_Num_Particles
use procedures_composition_factory, only: composition_destroy => destroy
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_different
use types_dipolar_interactions_dynamic_wrapper, only: Dipolar_Interactions_Dynamic_Wrapper
use procedures_dipoles_field_interaction, only: dipoles_field_visit_add => visit_add
use classes_random_coordinates, only: Abstract_Random_Coordinates
use procedures_random_coordinates_factory, only: random_coordinates_destroy => destroy
use types_observables_energies, only: Concrete_Single_Energies
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: &
        Abstract_Particle_Insertion_Method
    private
        type(Environment_Wrapper), pointer :: environment => null()
        class(Abstract_Num_Particles), allocatable :: nums_particles(:)
        type(Component_Wrapper), pointer :: components(:, :) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        type(Dipolar_Interactions_Dynamic_Wrapper), pointer :: dipolar_interactions_dynamic => &
            null()
        class(Abstract_Random_Coordinates), allocatable :: random_positions(:), random_orientation
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
        procedure, private :: visit_walls => Abstract_visit_walls
        procedure, private :: visit_short => Abstract_visit_short
        procedure, private :: visit_field => Abstract_visit_field
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

    subroutine Abstract_construct(this, environment, nums_particles, components, &
        short_interactions, dipolar_interactions_dynamic, random_positions, random_orientation)
        class(Abstract_Particle_Insertion_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Num_Particles), intent(in) :: nums_particles(:)
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        class(Abstract_Random_Coordinates), intent(in) :: random_positions(:), random_orientation

        this%environment => environment
        allocate(this%nums_particles, source=nums_particles)
        this%components => components
        this%short_interactions => short_interactions
        this%dipolar_interactions_dynamic => dipolar_interactions_dynamic
        allocate(this%random_positions, source=random_positions)
        allocate(this%random_orientation, source=random_orientation)
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Particle_Insertion_Method), intent(inout) :: this

        call random_coordinates_destroy(this%random_orientation)
        call random_coordinates_destroy(this%random_positions)
        call composition_destroy(this%nums_particles)
        this%dipolar_interactions_dynamic => null()
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
        integer :: i_box, i_component, i_particle
        logical :: overlap

        allocate(deltas%short_energies(size(observables%inv_pow_activities)))
        allocate(deltas%dipolar_energies(size(deltas%short_energies)))
        observables%gemc_inv_pow_activities = 0._DP
        test%i = 0
        do i_box = 1, size(this%environment%periodic_boxes)
            do i_component = 1, size(this%nums_particles)
                inv_pow_activity_sum = 0._DP
                do i_particle = 1, this%nums_particles(i_component)%get()

                    observables%gemc_insertion_counters(i_component, i_box)%num_hits = observables%&
                        gemc_insertion_counters(i_component, i_box)%num_hits + 1
                    test%position = this%random_positions(i_box)%get(i_component)
                    test%orientation = this%random_orientation%get(i_component)
                    test%dipole_moment = this%components(i_component, i_box)%dipole_moments%get_norm() * &
                        test%orientation

                    call this%visit_walls(overlap, deltas%walls_energy, i_box, i_component, test)
                    if (overlap) cycle
                    call this%visit_short(overlap, deltas%short_energies, i_box, i_component, test)
                    if (overlap) cycle
                    call this%visit_field(deltas%field_energy, i_box, test)
                    call this%visit_dipolar(deltas%dipolar_energies, deltas%dipolar_shared_energy, i_box, &
                        i_component, test)

                    delta_energy = deltas%field_energy + deltas%walls_energy + &
                        sum(deltas%short_energies + deltas%dipolar_energies) + deltas%&
                        dipolar_shared_energy
                    inv_pow_activity_sum = inv_pow_activity_sum + exp(-delta_energy / this%environment%&
                        temperature%get())
                    observables%gemc_insertion_counters(i_component, i_box)%num_successes = observables%&
                        gemc_insertion_counters(i_component, i_box)%num_successes + 1
                end do
                observables%gemc_inv_pow_activities(i_component, i_box) = inv_pow_activity_sum / real(this%&
                    nums_particles(i_component)%get())
            end do
        end do
    end subroutine Abstract_try

    subroutine Abstract_visit_walls(this, overlap, delta, i_box, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: delta
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        call this%environment%gemc_visitable_walls(i_box)%visit(overlap, delta, test%position, this%&
            short_interactions%wall_pairs(i_component)%potential)
    end subroutine Abstract_visit_walls

    subroutine Abstract_visit_short(this, overlap, deltas, i_box, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: deltas(:)
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        integer :: j_component, i_exclude

        do j_component = 1, size(this%short_interactions%cells(i_box)%visitable_cells, 1)
            i_exclude = merge(test%i, 0, j_component == i_component)
            call this%short_interactions%cells(i_box)%visitable_cells(j_component, i_component)%&
                visit_energy(overlap, deltas(j_component), test, visit_different, i_exclude)
                if (overlap) return
        end do
    end subroutine Abstract_visit_short

    subroutine Abstract_visit_field(this, delta, i_box, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        real(DP), intent(out) :: delta
        integer, intent(in) :: i_box
        type(Concrete_Temporary_Particle), intent(in) :: test

        delta = dipoles_field_visit_add(this%environment%external_fields(i_box), test)
    end subroutine Abstract_visit_field

    subroutine Abstract_visit_dipolar(this, deltas, delta_shared_energy, i_box, i_component, test)
        class(Abstract_Particle_Insertion_Method), intent(in) :: this
        real(DP), intent(out) :: deltas(:)
        real(DP), intent(out) :: delta_shared_energy
        integer, intent(in) :: i_box, i_component
        type(Concrete_Temporary_Particle), intent(in) :: test

        integer :: j_component, i_exclude

        do j_component = 1, size(this%dipolar_interactions_dynamic%gemc_real_components, 1)
            i_exclude = merge(test%i, 0, j_component == i_component)
            call this%dipolar_interactions_dynamic%gemc_real_components(j_component, i_component, i_box)%&
                component%visit(deltas(j_component), test, visit_different, i_exclude)
        end do
        deltas(i_component) = deltas(i_component) - this%dipolar_interactions_dynamic%&
            gemc_self_components(i_component, i_box)%component%meet(test%dipole_moment)
        delta_shared_energy = &
            this%dipolar_interactions_dynamic%reci_visitors(i_box)%visit_add(i_component, test) + &
            this%dipolar_interactions_dynamic%gemc_surf_mixture(i_box)%&
                visit_add(i_component, test%dipole_moment) - &
            this%dipolar_interactions_dynamic%dlc_visitors(i_box)%visit_add(i_component, test)
    end subroutine Abstract_visit_dipolar

!end implementation Abstract_Particle_Insertion_Method

!implementation Null_Particle_Insertion_Method

    subroutine Null_construct(this, environment, nums_particles, components, short_interactions, &
        dipolar_interactions_dynamic, random_positions, random_orientation)
        class(Null_Particle_Insertion_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        class(Abstract_Num_Particles), intent(in) :: nums_particles(:)
        type(Component_Wrapper), target, intent(in) :: components(:, :)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        type(Dipolar_Interactions_Dynamic_Wrapper), target, intent(in) :: &
            dipolar_interactions_dynamic
        class(Abstract_Random_Coordinates), intent(in) :: random_positions(:), random_orientation
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
