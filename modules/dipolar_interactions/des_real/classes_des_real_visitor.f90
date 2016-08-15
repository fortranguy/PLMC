module classes_des_real_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_component_coordinates, only: Abstract_Component_Coordinates
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_des_real_pair, only: Abstract_DES_Real_Pair

implicit none

private

    type, abstract, public :: Abstract_DES_Real_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        generic :: visit => visit_intra, visit_inter
        procedure, private :: visit_intra => Abstract_visit_intra
        procedure, private :: visit_inter => Abstract_visit_inter
    end type Abstract_DES_Real_Visitor

    type, extends(Abstract_DES_Real_Visitor), public :: Concrete_DES_Real_Visitor

    end type Concrete_DES_Real_Visitor

    type, extends(Abstract_DES_Real_Visitor), public :: Null_DES_Real_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure, private :: visit_intra => Null_visit_intra
        procedure, private :: visit_inter => Null_visit_inter
    end type Null_DES_Real_Visitor

contains

!implementation Abstract_DES_Real_Visitor

    subroutine Abstract_construct(this, periodic_box)
        class(Abstract_DES_Real_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Real_Visitor), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure subroutine Abstract_visit_inter(this, energy, positions_1, dipole_moments_1, positions_2, &
        dipole_moments_2, des_real_pair)
        class(Abstract_DES_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments_1, dipole_moments_2
        class(Abstract_DES_Real_Pair), intent(in) :: des_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: i_particle, j_particle

        energy = 0._DP
        if (dipole_moments_1%get_num() == 0 .or. dipole_moments_2%get_num() == 0) return
        do j_particle = 1, dipole_moments_2%get_num()
            do i_particle = 1, dipole_moments_1%get_num()
                vector_ij = this%periodic_box%vector(positions_1%get(i_particle), positions_2%&
                    get(j_particle))
                energy = energy + des_real_pair%meet(vector_ij, dipole_moments_1%&
                    get(i_particle), dipole_moments_2%get(j_particle))
            end do
        end do
    end subroutine Abstract_visit_inter

    pure subroutine Abstract_visit_intra(this, energy, positions, dipole_moments, des_real_pair)
        class(Abstract_DES_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), intent(in) :: des_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: i_particle, j_particle

        energy = 0._DP
        do j_particle = 1, dipole_moments%get_num()
            do i_particle = j_particle + 1, dipole_moments%get_num()
                vector_ij = this%periodic_box%vector(positions%get(i_particle), positions%&
                    get(j_particle))
                energy = energy + des_real_pair%meet(vector_ij, dipole_moments%get(i_particle), &
                    dipole_moments%get(j_particle))
            end do
        end do
    end subroutine Abstract_visit_intra

!end implementation Abstract_DES_Real_Visitor

!implementation Null_DES_Real_Visitor

    subroutine Null_construct(this, periodic_box)
        class(Null_DES_Real_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Real_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure subroutine Null_visit_inter(this, energy, positions_1, dipole_moments_1, positions_2, &
        dipole_moments_2, des_real_pair)
        class(Null_DES_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments_1, dipole_moments_2
        class(Abstract_DES_Real_Pair), intent(in) :: des_real_pair
        energy = 0._DP
    end subroutine Null_visit_inter

    pure subroutine Null_visit_intra(this, energy, positions, dipole_moments, des_real_pair)
        class(Null_DES_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipole_Moments), intent(in) :: dipole_moments
        class(Abstract_DES_Real_Pair), intent(in) :: des_real_pair
        energy = 0._DP
    end subroutine Null_visit_intra

!end implementation Null_DES_Real_Visitor

end module classes_des_real_visitor
