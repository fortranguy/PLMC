module classes_des_real_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle

implicit none

private

    type, abstract, public :: Abstract_DES_Real_Component
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Component_Dipole_Moments), pointer :: dipole_moments => null()
        class(Abstract_DES_Real_Pair), pointer :: des_real_pair => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: visit => visit_energy, visit_field
        procedure, private :: visit_energy => Abstract_visit_energy
        procedure, private :: visit_field => Abstract_visit_field
    end type Abstract_DES_Real_Component

    type, extends(Abstract_DES_Real_Component), public :: Concrete_DES_Real_Component

    end type Concrete_DES_Real_Component

    type, extends(Abstract_DES_Real_Component), public :: Null_DES_Real_Component
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure, private :: visit_energy => Null_visit_energy
        procedure, private :: visit_field => Null_visit_field
    end type Null_DES_Real_Component

contains

!implementation Abstract_DES_Real_Component

    subroutine Abstract_construct(this, periodic_box, positions, dipole_moments, des_real_pair)
        class(Abstract_DES_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), target, intent(in) :: des_real_pair

        this%periodic_box => periodic_box
        this%positions => positions
        this%dipole_moments => dipole_moments
        this%des_real_pair => des_real_pair
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Real_Component), intent(inout) :: this

        this%des_real_pair => null()
        this%dipole_moments => null()
        this%positions => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_visit_energy(this, energy, particle, i_exclude)
        class(Abstract_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        energy = 0._DP
        do j_particle = 1, this%positions%get_num()
            if (j_particle == i_exclude) cycle
            vector_ij = this%periodic_box%vector(particle%position, this%positions%&
                get(j_particle))
            energy = energy + this%des_real_pair%meet(vector_ij, particle%dipole_moment, &
                this%dipole_moments%get(j_particle))
        end do
    end subroutine Abstract_visit_energy

    pure subroutine Abstract_visit_field(this, field, particle, i_exclude)
        class(Abstract_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude

        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        field = 0._DP
        do j_particle = 1, this%positions%get_num()
            if (j_particle == i_exclude) cycle
            vector_ij = this%periodic_box%vector(particle%position, this%positions%get(j_particle))
            field = field + this%des_real_pair%meet(vector_ij, this%dipole_moments%get(j_particle))
        end do
    end subroutine Abstract_visit_field

!end implementation Abstract_DES_Real_Component

!implementation Null_DES_Real_Component

    subroutine Null_construct(this, periodic_box, positions, dipole_moments, des_real_pair)
        class(Null_DES_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), target, intent(in) :: des_real_pair
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Real_Component), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit_energy(this, energy, particle, i_exclude)
        class(Null_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude
        energy = 0._DP
    end subroutine Null_visit_energy

    pure subroutine Null_visit_field(this, field, particle, i_exclude)
        class(Null_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: field(num_dimensions)
        type(Concrete_Temporary_Particle), intent(in) :: particle
        integer, intent(in) :: i_exclude
        field = 0._DP
    end subroutine Null_visit_field

!end implementation Null_DES_Real_Component

end module classes_des_real_component
