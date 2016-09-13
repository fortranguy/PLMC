module classes_des_self_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use classes_permittivity, only: Abstract_Permittivity
use classes_component_dipole_moments, only: Abstract_Component_Dipole_Moments
use classes_des_convergence_parameter, only: Abstract_DES_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_DES_Self_Component
    private
        real(DP) :: permittivity = 0._DP
        class(Abstract_DES_Convergence_Parameter), pointer :: alpha => null()
        real(DP) :: alpha_cube = 0._DP
        real(DP) :: factor_denominator = 1._DP
        class(Abstract_Component_Dipole_Moments), pointer :: dipole_moments => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_set
        procedure :: visit => Abstract_visit
        procedure :: meet => Abstract_meet
        procedure, private :: set => Abstract_set
    end type Abstract_DES_Self_Component

    type, extends(Abstract_DES_Self_Component), public :: Concrete_DES_Self_Component

    end type Concrete_DES_Self_Component

    type, extends(Abstract_DES_Self_Component), public :: Null_DES_Self_Component
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_set
        procedure :: visit => Null_visit
        procedure :: meet => Null_meet
    end type Null_DES_Self_Component

contains

!implementation Abstract_DES_Self_Component

    subroutine Abstract_construct(this, permittivity, dipole_moments, alpha)
        class(Abstract_DES_Self_Component), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha

        this%permittivity = permittivity%get()
        this%dipole_moments => dipole_moments
        this%alpha => alpha
        call this%set()
    end subroutine Abstract_construct

    subroutine Abstract_set(this)
        class(Abstract_DES_Self_Component), intent(inout) :: this

        this%alpha_cube = this%alpha%get()**3
        this%factor_denominator = 6._DP * this%permittivity * PI**(3._DP/2._DP)
    end subroutine Abstract_set

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Self_Component), intent(inout) :: this

        this%alpha => null()
        this%dipole_moments => null()
    end subroutine Abstract_destroy

    !> \[ \sum_{\vec{\mu}} u(\vec{\mu}) \]
    !> where \( u(\vec{\mu}) \) is [[Abstract_meet]].
    pure real(DP) function Abstract_visit(this) result(energy)
        class(Abstract_DES_Self_Component), intent(in) :: this

        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, this%dipole_moments%get_num()
            energy = energy + this%meet(this%dipole_moments%get(i_particle))
        end do
    end function Abstract_visit

    !> \[ u(\mu) = \frac{\alpha^3}{6\epsilon\pi^{3/2}} \vec{\mu}\cdot\vec{\mu} \]
    pure real(DP) function Abstract_meet(this, dipole_moment) result(energy)
        class(Abstract_DES_Self_Component), intent(in) :: this
        real(DP), intent(in) :: dipole_moment(:)

        energy = this%alpha_cube/this%factor_denominator * dot_product(dipole_moment, dipole_moment)
    end function Abstract_meet

!end implementation Abstract_DES_Self_Component

!implementation Null_DES_Self_Component

    subroutine Null_construct(this, permittivity, dipole_moments, alpha)
        class(Null_DES_Self_Component), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Component_Dipole_Moments), target, intent(in) :: dipole_moments
        class(Abstract_DES_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_construct

    subroutine Null_set(this)
        class(Null_DES_Self_Component), intent(inout) :: this
    end subroutine Null_set

    subroutine Null_destroy(this)
        class(Null_DES_Self_Component), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_DES_Self_Component), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_meet(this, dipole_moment) result(self)
        class(Null_DES_Self_Component), intent(in) :: this
        real(DP), intent(in) :: dipole_moment(:)
        self = 0._DP
    end function Null_meet

!end implementation Null_DES_Self_Component

end module classes_des_self_component
