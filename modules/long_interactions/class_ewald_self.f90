module class_ewald_self

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: PI
use class_permittivity, only: Abstract_Permittivity
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter

implicit none

private

    type, abstract, public :: Abstract_Ewald_Self
    private
        class(Abstract_Permittivity), pointer :: permittivity => null()
        class(Abstract_Ewald_Convergence_Parameter), pointer :: alpha => null()
        real(DP) :: apha_cube
        real(DP) :: factor_denominator
    contains
        procedure :: construct => Abstract_Ewald_Self_construct
        procedure :: destroy => Abstract_Ewald_Self_destroy
        procedure :: reset => Abstract_Ewald_Self_set
        procedure :: get => Abstract_Ewald_Self_get
        procedure, private :: set => Abstract_Ewald_Self_set
    end type Abstract_Ewald_Self

    type, extends(Abstract_Ewald_Self), public :: Concrete_Ewald_Self

    end type Concrete_Ewald_Self

    type, extends(Abstract_Ewald_Self), public :: Null_Ewald_Self
    contains
        procedure :: construct => Null_Ewald_Self_construct
        procedure :: destroy => Null_Ewald_Self_destroy
        procedure :: reset => Null_Ewald_Self_set
        procedure :: get => Null_Ewald_Self_get
    end type Null_Ewald_Self

contains

!implementation Abstract_Ewald_Self

    subroutine Abstract_Ewald_Self_construct(this, permittivity, alpha)
        class(Abstract_Ewald_Self), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha

        this%permittivity => permittivity
        this%alpha => alpha
        call this%set()
    end subroutine Abstract_Ewald_Self_construct

    subroutine Abstract_Ewald_Self_set(this)
        class(Abstract_Ewald_Self), intent(inout) :: this

        this%apha_cube = this%alpha%get()**3
        this%factor_denominator = 6._DP * this%permittivity%get() * PI**(3._DP/2._DP)
    end subroutine Abstract_Ewald_Self_set

    subroutine Abstract_Ewald_Self_destroy(this)
        class(Abstract_Ewald_Self), intent(inout) :: this

        this%alpha => null()
        this%permittivity => null()
    end subroutine Abstract_Ewald_Self_destroy

    !> \[ u(\mu) = \frac{\alpha^3}{6\epsilon\pi^{3/2}} \vec{\mu}\cdot\vec{\mu} \]
    pure real(DP) function Abstract_Ewald_Self_get(this, moment) result(self)
        class(Abstract_Ewald_Self), intent(in) :: this
        real(DP), intent(in) :: moment(:)

        self = this%apha_cube/this%factor_denominator * dot_product(moment, moment)
    end function Abstract_Ewald_Self_get

!end implementation Abstract_Ewald_Self

!implementation Null_Ewald_Self

    subroutine Null_Ewald_Self_construct(this, permittivity, alpha)
        class(Null_Ewald_Self), intent(out) :: this
        class(Abstract_Permittivity), target, intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), target, intent(in) :: alpha
    end subroutine Null_Ewald_Self_construct

    subroutine Null_Ewald_Self_set(this)
        class(Null_Ewald_Self), intent(inout) :: this
    end subroutine Null_Ewald_Self_set

    subroutine Null_Ewald_Self_destroy(this)
        class(Null_Ewald_Self), intent(inout) :: this
    end subroutine Null_Ewald_Self_destroy

    pure real(DP) function Null_Ewald_Self_get(this, moment) result(self)
        class(Null_Ewald_Self), intent(in) :: this
        real(DP), intent(in) :: moment(:)
        self = 0._DP
    end function Null_Ewald_Self_get

!end implementation Null_Ewald_Self

end module class_ewald_self
