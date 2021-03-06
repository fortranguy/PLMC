module classes_short_pairs_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_pair_potential, only: Abstract_Pair_Potential

implicit none

private

    type, abstract, public :: Abstract_Short_Pairs_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: visit => visit_intra, visit_inter
        procedure, private :: visit_intra => Abstract_visit_intra
        procedure, private :: visit_inter => Abstract_visit_inter
    end type Abstract_Short_Pairs_Visitor

    type, extends(Abstract_Short_Pairs_Visitor), public :: Concrete_Short_Pairs_Visitor

    end type Concrete_Short_Pairs_Visitor

    type, extends(Abstract_Short_Pairs_Visitor), public :: Null_Short_Pairs_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure, private :: visit_intra => Null_visit_intra
        procedure, private :: visit_inter => Null_visit_inter
    end type Null_Short_Pairs_Visitor

contains

!implementation Abstract_Short_Pairs_Visitor

    subroutine Abstract_construct(this, periodic_box)
        class(Abstract_Short_Pairs_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Short_Pairs_Visitor), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_visit_inter(this, overlap, energy, positions_1, positions_2, &
        pair_potential)
        class(Abstract_Short_Pairs_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_ij, distance_ij
        integer :: i_particle, j_particle

        overlap = .false.
        energy = 0._DP
        do j_particle = 1, positions_2%get_num()
            do i_particle = 1, positions_1%get_num()
                distance_ij = this%periodic_box%distance(positions_1%get(i_particle), positions_2%&
                    get(j_particle))
                call pair_potential%meet(overlap, energy_ij, distance_ij)
                if (overlap) return
                energy = energy + energy_ij
            end do
        end do
    end subroutine Abstract_visit_inter

    pure subroutine Abstract_visit_intra(this, overlap, energy, positions, pair_potential)
        class(Abstract_Short_Pairs_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        real(DP) :: energy_ij, distance_ij
        integer :: i_particle, j_particle

        overlap = .false.
        energy = 0._DP
        do j_particle = 1, positions%get_num()
            do i_particle = j_particle + 1, positions%get_num()
                distance_ij = this%periodic_box%distance(positions%get(i_particle), positions%&
                    get(j_particle))
                call pair_potential%meet(overlap, energy_ij, distance_ij)
                if (overlap) return
                energy = energy + energy_ij
            end do
        end do
    end subroutine Abstract_visit_intra

!end implementation Abstract_Short_Pairs_Visitor

!implementation Null_Short_Pairs_Visitor

    subroutine Null_construct(this, periodic_box)
        class(Null_Short_Pairs_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Short_Pairs_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit_inter(this, overlap, energy, positions_1, positions_2, &
        pair_potential)
        class(Null_Short_Pairs_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit_inter

    pure subroutine Null_visit_intra(this, overlap, energy, positions, pair_potential)
        class(Null_Short_Pairs_Visitor), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Pair_Potential), intent(in) :: pair_potential
        overlap = .false.
        energy = 0._DP
    end subroutine Null_visit_intra

!end implementation Null_Short_Pairs_Visitor

end module classes_short_pairs_visitor
