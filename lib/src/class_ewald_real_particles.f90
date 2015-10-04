module class_ewald_real_particles

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_periodic_box, only: Abstract_Periodic_Box
use class_particles_positions, only: Abstract_Particles_Positions
use class_particles_dipolar_moments, only: Abstract_Particles_Dipolar_Moments
use types_particle, only: Concrete_Particle

implicit none

private

    type, abstract, public :: Abstract_Ewald_Real_Particles
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Particles_Positions), pointer :: particles_positions => null()
        class(Abstract_Particles_Dipolar_Moments), pointer :: particles_dipolar_moments => null()
    contains
        procedure :: construct => Abstract_Ewald_Real_Particles_construct
        procedure :: destroy => Abstract_Ewald_Real_Particles_destroy
        generic :: visit => visit_energy, visit_field
        procedure, private :: visit_energy => Abstract_Ewald_Real_Particles_visit_energy
        procedure, private :: visit_field => Abstract_Ewald_Real_Particles_visit_field
    end type Abstract_Ewald_Real_Particles

    type, extends(Abstract_Ewald_Real_Particles) :: Concrete_Ewald_Real_Particles

    end type Concrete_Ewald_Real_Particles

    type, extends(Abstract_Ewald_Real_Particles), public :: Null_Ewald_Real_Particles
    contains
        procedure :: construct => Null_Ewald_Real_Particles_construct
        procedure :: destroy => Null_Ewald_Real_Particles_destroy
        procedure, private :: visit_energy => Null_Ewald_Real_Particles_visit_energy
        procedure, private :: visit_field => Null_Ewald_Real_Particles_visit_field
    end type Null_Ewald_Real_Particles

contains

!implementation Abstract_Ewald_Real_Particles

    subroutine Abstract_Ewald_Real_Particles_construct(this, periodic_box, particles_positions, &
        particles_dipolar_moments)
        class(Abstract_Ewald_Real_Particles), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Positions), target, intent(in) :: particles_positions
        class(Abstract_Particles_Dipolar_Moments), target, intent(in) :: particles_dipolar_moments

        this%periodic_box => periodic_box
        this%particles_positions => particles_positions
        this%particles_dipolar_moments => particles_dipolar_moments
    end subroutine Abstract_Ewald_Real_Particles_construct

    subroutine Abstract_Ewald_Real_Particles_destroy(this)
        class(Abstract_Ewald_Real_Particles), intent(inout) :: this

        this%particles_dipolar_moments => null()
        this%particles_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Real_Particles_destroy

    pure subroutine Abstract_Ewald_Real_Particles_visit_energy(this, energy, particle, &
        ewald_real_pair)
        class(Abstract_Ewald_Real_Particles), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        energy = 0._DP
        do j_particle = 1, this%particles_positions%get_num()
            if (particle%same_type .and. particle%i == j_particle) cycle
            vector_ij = this%periodic_box%vector(particle%position, &
                this%particles_positions%get(j_particle))
            energy = energy + ewald_real_pair%meet(vector_ij, particle%dipolar_moment, &
                this%particles_dipolar_moments%get(j_particle))
        end do
    end subroutine Abstract_Ewald_Real_Particles_visit_energy

    pure subroutine Abstract_Ewald_Real_Particles_visit_field(this, field, particle, &
        ewald_real_pair)
        class(Abstract_Ewald_Real_Particles), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        field = 0._DP
        do j_particle = 1, this%particles_positions%get_num()
            if (particle%same_type .and. particle%i == j_particle) cycle
            vector_ij = this%periodic_box%vector(particle%position, &
                this%particles_positions%get(j_particle))
            field = field + ewald_real_pair%meet(vector_ij, &
                this%particles_dipolar_moments%get(j_particle))
        end do
    end subroutine Abstract_Ewald_Real_Particles_visit_field

!end implementation Abstract_Ewald_Real_Particles

!implementation Null_Ewald_Real_Particles

    subroutine Null_Ewald_Real_Particles_construct(this, periodic_box, particles_positions, &
        particles_dipolar_moments)
        class(Null_Ewald_Real_Particles), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Particles_Positions), target, intent(in) :: particles_positions
        class(Abstract_Particles_Dipolar_Moments), target, intent(in) :: particles_dipolar_moments
    end subroutine Null_Ewald_Real_Particles_construct

    subroutine Null_Ewald_Real_Particles_destroy(this)
        class(Null_Ewald_Real_Particles), intent(inout) :: this
    end subroutine Null_Ewald_Real_Particles_destroy

    pure subroutine Null_Ewald_Real_Particles_visit_energy(this, energy, particle, ewald_real_pair)
        class(Null_Ewald_Real_Particles), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        energy = 0._DP
    end subroutine Null_Ewald_Real_Particles_visit_energy

    pure subroutine Null_Ewald_Real_Particles_visit_field(this, field, particle, ewald_real_pair)
        class(Null_Ewald_Real_Particles), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        type(Concrete_Particle), intent(in) :: particle
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        field = 0._DP
    end subroutine Null_Ewald_Real_Particles_visit_field

!end implementation Null_Ewald_Real_Particles

end module class_ewald_real_particles
