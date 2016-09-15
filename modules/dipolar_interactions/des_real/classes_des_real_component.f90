module classes_des_real_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_des_real_pair, only: Abstract_DES_Real_Pair
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle
use procedures_visit_condition, only: abstract_visit_condition

implicit none

private

    type, abstract, public :: Abstract_DES_Real_Component
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Box_Volume_Memento), pointer :: box_volume_memento => null()
        class(Abstract_Component_Coordinates), pointer :: positions => null()
        class(Abstract_Component_Dipole_Moments), pointer :: dipole_moments => null()
        class(Abstract_DES_Real_Pair), pointer :: des_real_pair => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
    end type Abstract_DES_Real_Component

    type, extends(Abstract_DES_Real_Component), public :: Concrete_DES_Real_Component

    end type Concrete_DES_Real_Component

    type, extends(Abstract_DES_Real_Component), public :: Null_DES_Real_Component
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
    end type Null_DES_Real_Component

contains

!implementation Abstract_DES_Real_Component

    subroutine Abstract_construct(this, periodic_box, box_volume_memento, positions, &
        dipole_moments, des_real_pair)
        class(Abstract_DES_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Box_Volume_Memento), target, intent(in) :: box_volume_memento
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), target, intent(in) :: des_real_pair

        this%periodic_box => periodic_box
        this%box_volume_memento => box_volume_memento
        this%positions => positions
        this%dipole_moments => dipole_moments
        this%des_real_pair => des_real_pair
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Real_Component), intent(inout) :: this

        this%des_real_pair => null()
        this%dipole_moments => null()
        this%positions => null()
        this%box_volume_memento => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> \[
    !>      \frac{V_\text{s}}{V} \sum_{j} [\text{c}(j, i_\text{exclude})] u \left(\left(
    !>          \frac{V_\text{s}}{V} \right)^{1/3} \vec{r}_{ij}, \vec{\mu}_i, \vec{\mu}_j
    !>      \right)
    !> \]
    !> cf. [[classes_box_volume_memento:Abstract_get]] and
    !> [[classes_des_real_pair:Abstract_meet]]
    pure subroutine Abstract_visit(this, energy, particle, visit_condition, i_exclude)
        class(Abstract_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude

        real(DP) :: box_volume_ratio, box_edge_ratio
        real(DP) :: vector_ij(num_dimensions)
        integer :: j_particle

        box_volume_ratio = this%box_volume_memento%get() / product(this%periodic_box%get_size())
        box_edge_ratio = box_volume_ratio**(1._DP/3._DP)
        energy = 0._DP
        do j_particle = 1, this%positions%get_num()
            if (.not.visit_condition(j_particle, i_exclude)) cycle
            vector_ij = this%periodic_box%vector(particle%position, this%positions%get(j_particle))
            energy = energy + this%des_real_pair%meet(box_edge_ratio * vector_ij, particle%&
                dipole_moment, this%dipole_moments%get(j_particle))
        end do
        energy = box_volume_ratio * energy
    end subroutine Abstract_visit

!end implementation Abstract_DES_Real_Component

!implementation Null_DES_Real_Component

    subroutine Null_construct(this, periodic_box, box_volume_memento, positions, dipole_moments, &
        des_real_pair)
        class(Null_DES_Real_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Box_Volume_Memento), target, intent(in) :: box_volume_memento
        class(Abstract_Component_Coordinates), target, intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), target, intent(in) :: des_real_pair
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Real_Component), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit(this, energy, particle, visit_condition, i_exclude)
        class(Null_DES_Real_Component), intent(in) :: this
        real(DP), intent(out) :: energy
        type(Concrete_Temporary_Particle), intent(in) :: particle
        procedure(abstract_visit_condition) :: visit_condition
        integer, intent(in) :: i_exclude
        energy = 0._DP
    end subroutine Null_visit

!end implementation Null_DES_Real_Component

end module classes_des_real_component
