module class_ewald_real_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

    type, abstract, public :: Abstract_Ewald_Real_Component
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: component_positions => null()
        class(Abstract_Component_Dipolar_Moments), pointer :: component_dipolar_moments => null()
        class(Abstract_Ewald_Real_Pair), pointer :: ewald_real_pair => null()
    contains
        procedure :: construct => Abstract_Ewald_Real_Component_construct
        procedure :: destroy => Abstract_Ewald_Real_Component_destroy
        generic :: visit => visit_energy, visit_field
        procedure, private :: visit_energy => Abstract_Ewald_Real_Component_visit_energy
        procedure, private :: visit_field => Abstract_Ewald_Real_Component_visit_field
    end type Abstract_Ewald_Real_Component

    type, extends(Abstract_Ewald_Real_Component), public :: Concrete_Ewald_Real_Component

    end type Concrete_Ewald_Real_Component

    type, extends(Abstract_Ewald_Real_Component), public :: Null_Ewald_Real_Component
    contains
        procedure :: construct => Null_Ewald_Real_Component_construct
        procedure :: destroy => Null_Ewald_Real_Component_destroy
        procedure, private :: visit_energy => Null_Ewald_Real_Component_visit_energy
        procedure, private :: visit_field => Null_Ewald_Real_Component_visit_field
    end type Null_Ewald_Real_Component

contains

!implementation Abstract_Ewald_Real_Component

    subroutine Abstract_Ewald_Real_Component_construct(this, periodic_box, component_positions, &
        component_dipolar_moments)
        class(Abstract_Ewald_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments

        this%periodic_box => periodic_box
        this%component_positions => component_positions
        this%component_dipolar_moments => component_dipolar_moments
    end subroutine Abstract_Ewald_Real_Component_construct

    subroutine Abstract_Ewald_Real_Component_destroy(this)
        class(Abstract_Ewald_Real_Component), intent(inout) :: this

        this%component_dipolar_moments => null()
        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Real_Component_destroy

    pure subroutine Abstract_Ewald_Real_Component_visit_energy(this, energy, ewald_real_pair, &
        particle, i_exclude)
        class(Abstract_Ewald_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        energy = 0._DP
        do j_particle = 1, this%component_positions%get_num()
            if (j_particle == i_exclude) cycle
            vector_ij = this%periodic_box%vector(particle%position, this%component_positions%&
                get(j_particle))
            energy = energy + ewald_real_pair%meet(vector_ij, particle%dipolar_moment, &
                this%component_dipolar_moments%get(j_particle))
        end do
    end subroutine Abstract_Ewald_Real_Component_visit_energy

    pure subroutine Abstract_Ewald_Real_Component_visit_field(this, field, ewald_real_pair, &
        particle, i_exclude)
        class(Abstract_Ewald_Real_Component), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        field = 0._DP
        do j_particle = 1, this%component_positions%get_num()
            if (j_particle == i_exclude) cycle
            vector_ij = this%periodic_box%vector(particle%position, this%component_positions%&
                get(j_particle))
            field = field + ewald_real_pair%meet(vector_ij, this%component_dipolar_moments%&
                get(j_particle))
        end do
    end subroutine Abstract_Ewald_Real_Component_visit_field

!end implementation Abstract_Ewald_Real_Component

!implementation Null_Ewald_Real_Component

    subroutine Null_Ewald_Real_Component_construct(this, periodic_box, component_positions, &
        component_dipolar_moments)
        class(Null_Ewald_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments
    end subroutine Null_Ewald_Real_Component_construct

    subroutine Null_Ewald_Real_Component_destroy(this)
        class(Null_Ewald_Real_Component), intent(inout) :: this
    end subroutine Null_Ewald_Real_Component_destroy

    pure subroutine Null_Ewald_Real_Component_visit_energy(this, energy, ewald_real_pair, &
        particle, i_exclude)
        class(Null_Ewald_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude
        energy = 0._DP
    end subroutine Null_Ewald_Real_Component_visit_energy

    pure subroutine Null_Ewald_Real_Component_visit_field(this, field, ewald_real_pair, particle, &
        i_exclude)
        class(Null_Ewald_Real_Component), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude
        field = 0._DP
    end subroutine Null_Ewald_Real_Component_visit_field

!end implementation Null_Ewald_Real_Component

end module class_ewald_real_component
