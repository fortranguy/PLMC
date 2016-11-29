module classes_des_self_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Self_Component
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Component_Dipole_Moments), pointer :: dipole_moments => null()
        real(DP) :: factor_denominator = 0._DP
        real(DP) :: alpha_x_box_edge = 0._DP
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
        procedure :: meet => Abstract_meet
    end type Abstract_DES_Self_Component

    type, extends(Abstract_DES_Self_Component), public :: Concrete_DES_Self_Component

    end type Concrete_DES_Self_Component

    type, extends(Abstract_DES_Self_Component), public :: Null_DES_Self_Component
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
        procedure :: meet => Null_meet
    end type Null_DES_Self_Component

    type, public :: DES_Self_Component_Wrapper
        class(Abstract_DES_Self_Component), allocatable :: component
    end type DES_Self_Component_Wrapper

contains

!implementation Abstract_DES_Self_Component

    subroutine Abstract_construct(this, periodic_box, permittivity, dipole_moments, alpha)
        class(Abstract_DES_Self_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha

        this%periodic_box => periodic_box
        this%dipole_moments => dipole_moments
        this%factor_denominator = 6._DP * permittivity%get() * PI**(3._DP/2._DP)
        this%alpha_x_box_edge = alpha%get_times_box_edge()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Self_Component), intent(inout) :: this

        this%dipole_moments => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> \[ \frac{\alpha^3}{6\epsilon\pi^{3/2}} \sum_{\vec{\mu}} \vec{\mu}\cdot\vec{\mu} \]
    pure real(DP) function Abstract_visit(this) result(energy)
        class(Abstract_DES_Self_Component), intent(in) :: this

        real(DP) :: box_size(num_dimensions)
        real(DP) :: moment_squared_sum
        integer :: i_particle

        moment_squared_sum = 0._DP
        do i_particle = 1, this%dipole_moments%get_num()
            moment_squared_sum = moment_squared_sum + dot_product(this%dipole_moments%&
                get(i_particle), this%dipole_moments%get(i_particle))
        end do
        box_size = this%periodic_box%get_size()
        energy = (this%alpha_x_box_edge/box_size(1))**3 / this%factor_denominator * &
            moment_squared_sum
    end function Abstract_visit

    !> \[ u(\mu) = \frac{\alpha^3}{6\epsilon\pi^{3/2}} \vec{\mu}\cdot\vec{\mu} \]
    pure real(DP) function Abstract_meet(this, dipole_moment) result(energy)
        class(Abstract_DES_Self_Component), intent(in) :: this
        real(DP), intent(in) :: dipole_moment(:)

        real(DP) :: box_size(num_dimensions)

        box_size = this%periodic_box%get_size()
        energy = (this%alpha_x_box_edge/box_size(1))**3 / this%factor_denominator * &
            dot_product(dipole_moment, dipole_moment)
    end function Abstract_meet

!end implementation Abstract_DES_Self_Component

!implementation Null_DES_Self_Component

    subroutine Null_construct(this, periodic_box, permittivity, dipole_moments, alpha)
        class(Null_DES_Self_Component), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Convergence_Parameter), intent(in) :: alpha
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Self_Component), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_DES_Self_Component), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_meet(this, dipole_moment) result(energy)
        class(Null_DES_Self_Component), intent(in) :: this
        real(DP), intent(in) :: dipole_moment(:)
        energy = 0._DP
    end function Null_meet

!end implementation Null_DES_Self_Component

end module classes_des_self_component
