module class_ewald_real_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_ewald_real_pair, only: Abstract_Ewald_Real_Pair

implicit none

private

    type, abstract, public :: Abstract_Ewald_Real_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
    contains
        procedure :: construct => Abstract_Ewald_Real_Visitor_construct
        procedure :: destroy => Abstract_Ewald_Real_Visitor_destroy
        generic :: visit => visit_intra, visit_inter
        procedure, private :: visit_intra => Abstract_Ewald_Real_Visitor_visit_intra
        procedure, private :: visit_inter => Abstract_Ewald_Real_Visitor_visit_inter
    end type Abstract_Ewald_Real_Visitor

    type, extends(Abstract_Ewald_Real_Visitor), public :: Concrete_Ewald_Real_Visitor

    end type Concrete_Ewald_Real_Visitor

    type, extends(Abstract_Ewald_Real_Visitor), public :: Null_Ewald_Real_Visitor
    contains
        procedure :: construct => Null_Ewald_Real_Visitor_construct
        procedure :: destroy => Null_Ewald_Real_Visitor_destroy
        procedure, private :: visit_intra => Null_Ewald_Real_Visitor_visit_intra
        procedure, private :: visit_inter => Null_Ewald_Real_Visitor_visit_inter
    end type Null_Ewald_Real_Visitor

contains

!implementation Abstract_Ewald_Real_Visitor

    subroutine Abstract_Ewald_Real_Visitor_construct(this, periodic_box)
        class(Abstract_Ewald_Real_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box

        this%periodic_box => periodic_box
    end subroutine Abstract_Ewald_Real_Visitor_construct

    subroutine Abstract_Ewald_Real_Visitor_destroy(this)
        class(Abstract_Ewald_Real_Visitor), intent(inout) :: this

        this%periodic_box => null()
    end subroutine Abstract_Ewald_Real_Visitor_destroy

    pure subroutine Abstract_Ewald_Real_Visitor_visit_inter(this, energy, positions_1, &
        dipolar_moment_1, positions_2, dipolar_moment_2, ewald_real_pair)
        class(Abstract_Ewald_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moment_1, dipolar_moment_2
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: i_particle, j_particle

        energy = 0._DP
        do j_particle = 1, dipolar_moment_2%get_num()
            do i_particle = 1, dipolar_moment_1%get_num()
                vector_ij = this%periodic_box%vector(positions_1%get(i_particle), positions_2%&
                    get(j_particle))
                energy = energy + ewald_real_pair%meet(vector_ij, dipolar_moment_1%&
                    get(i_particle), dipolar_moment_2%get(j_particle))
            end do
        end do
    end subroutine Abstract_Ewald_Real_Visitor_visit_inter

    pure subroutine Abstract_Ewald_Real_Visitor_visit_intra(this, energy, positions, &
        dipolar_moment, ewald_real_pair)
        class(Abstract_Ewald_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moment
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair

        real(DP) :: vector_ij(num_dimensions)
        integer :: i_particle, j_particle

        energy = 0._DP
        do j_particle = 1, dipolar_moment%get_num()
            do i_particle = j_particle + 1, dipolar_moment%get_num()
                vector_ij = this%periodic_box%vector(positions%get(i_particle), positions%&
                    get(j_particle))
                energy = energy + ewald_real_pair%meet(vector_ij, dipolar_moment%get(i_particle), &
                    dipolar_moment%get(j_particle))
            end do
        end do
    end subroutine Abstract_Ewald_Real_Visitor_visit_intra

!end implementation Abstract_Ewald_Real_Visitor

!implementation Null_Ewald_Real_Visitor

    subroutine Null_Ewald_Real_Visitor_construct(this, periodic_box)
        class(Null_Ewald_Real_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
    end subroutine Null_Ewald_Real_Visitor_construct

    subroutine Null_Ewald_Real_Visitor_destroy(this)
        class(Null_Ewald_Real_Visitor), intent(inout) :: this
    end subroutine Null_Ewald_Real_Visitor_destroy

    pure subroutine Null_Ewald_Real_Visitor_visit_inter(this, energy, positions_1, &
        dipolar_moment_1, positions_2, dipolar_moment_2, ewald_real_pair)
        class(Null_Ewald_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions_1, positions_2
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moment_1, dipolar_moment_2
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        energy = 0._DP
    end subroutine Null_Ewald_Real_Visitor_visit_inter

    pure subroutine Null_Ewald_Real_Visitor_visit_intra(this, energy, positions, dipolar_moment, &
        ewald_real_pair)
        class(Null_Ewald_Real_Visitor), intent(in) :: this
        real(DP), intent(out) :: energy
        class(Abstract_Component_Coordinates), intent(in) :: positions
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moment
        class(Abstract_Ewald_Real_Pair), intent(in) :: ewald_real_pair
        energy = 0._DP
    end subroutine Null_Ewald_Real_Visitor_visit_intra

!end implementation Null_Ewald_Real_Visitor

end module class_ewald_real_visitor
