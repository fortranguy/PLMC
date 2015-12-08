module class_ewald_self_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use class_permittivity, only: Abstract_Permittivity
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_Ewald_Self_Component
    private
        class(Abstract_Permittivity), pointer :: permittivity => null()
        class(Abstract_Ewald_Convergence_Parameter), pointer :: alpha => null()
        real(DP) :: apha_cube
        real(DP) :: factor_denominator
        class(Abstract_Component_Dipolar_Moments), pointer :: component_dipolar_moment => null()
    contains
        procedure :: construct => Abstract_Ewald_Self_Component_construct
        procedure :: destroy => Abstract_Ewald_Self_Component_destroy
        procedure :: reset => Abstract_Ewald_Self_Component_set
        procedure :: meet => Abstract_Ewald_Self_Component_meet
        procedure :: visit => Abstract_Ewald_Self_Component_visit
        procedure, private :: set => Abstract_Ewald_Self_Component_set
    end type Abstract_Ewald_Self_Component

    type, extends(Abstract_Ewald_Self_Component), public :: Concrete_Ewald_Self_Component

    end type Concrete_Ewald_Self_Component

    type, extends(Abstract_Ewald_Self_Component), public :: Null_Ewald_Self_Component
    contains
        procedure :: construct => Null_Ewald_Self_Component_construct
        procedure :: destroy => Null_Ewald_Self_Component_destroy
        procedure :: reset => Null_Ewald_Self_Component_set
        procedure :: visit => Null_Ewald_Self_Component_visit
        procedure :: meet => Null_Ewald_Self_Component_meet
    end type Null_Ewald_Self_Component

contains

!implementation Abstract_Ewald_Self_Component

    subroutine Abstract_Ewald_Self_Component_construct(this, permittivity, &
        component_dipolar_moment, alpha)
        class(Abstract_Ewald_Self_Component), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moment
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha

        this%permittivity => permittivity
        this%component_dipolar_moment => component_dipolar_moment
        this%alpha => alpha
        call this%set()
    end subroutine Abstract_Ewald_Self_Component_construct

    subroutine Abstract_Ewald_Self_Component_set(this)
        class(Abstract_Ewald_Self_Component), intent(inout) :: this

        this%apha_cube = this%alpha%get()**3
        this%factor_denominator = 6._DP * this%permittivity%get() * PI**(3._DP/2._DP)
    end subroutine Abstract_Ewald_Self_Component_set

    subroutine Abstract_Ewald_Self_Component_destroy(this)
        class(Abstract_Ewald_Self_Component), intent(inout) :: this

        this%alpha => null()
        this%component_dipolar_moment => null()
        this%permittivity => null()
    end subroutine Abstract_Ewald_Self_Component_destroy

    pure real(DP) function Abstract_Ewald_Self_Component_visit(this) result(energy)
        class(Abstract_Ewald_Self_Component), intent(in) :: this

        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, this%component_dipolar_moment%get_num()
            energy = energy + this%meet(this%component_dipolar_moment%get(i_particle))
        end do
    end function Abstract_Ewald_Self_Component_visit

    !> \[ u(\mu) = \frac{\alpha^3}{6\epsilon\pi^{3/2}} \vec{\mu}\cdot\vec{\mu} \]
    pure real(DP) function Abstract_Ewald_Self_Component_meet(this, moment) result(energy)
        class(Abstract_Ewald_Self_Component), intent(in) :: this
        real(DP), intent(in) :: moment(:)

        energy = this%apha_cube/this%factor_denominator * dot_product(moment, moment)
    end function Abstract_Ewald_Self_Component_meet

!end implementation Abstract_Ewald_Self_Component

!implementation Null_Ewald_Self_Component

    subroutine Null_Ewald_Self_Component_construct(this, permittivity, component_dipolar_moment, &
        alpha)
        class(Null_Ewald_Self_Component), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moment
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_Ewald_Self_Component_construct

    subroutine Null_Ewald_Self_Component_set(this)
        class(Null_Ewald_Self_Component), intent(inout) :: this
    end subroutine Null_Ewald_Self_Component_set

    subroutine Null_Ewald_Self_Component_destroy(this)
        class(Null_Ewald_Self_Component), intent(inout) :: this
    end subroutine Null_Ewald_Self_Component_destroy

    pure real(DP) function Null_Ewald_Self_Component_visit(this) result(energy)
        class(Null_Ewald_Self_Component), intent(in) :: this
        energy = 0._DP
    end function Null_Ewald_Self_Component_visit

    pure real(DP) function Null_Ewald_Self_Component_meet(this, moment) result(self)
        class(Null_Ewald_Self_Component), intent(in) :: this
        real(DP), intent(in) :: moment(:)
        self = 0._DP
    end function Null_Ewald_Self_Component_meet

!end implementation Null_Ewald_Self_Component

end module class_ewald_self_component
